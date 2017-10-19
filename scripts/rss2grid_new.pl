#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module
use Statistics::OLS;
use Math::FFT;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);
use Math::Approx;
use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Slatec;
use PDL::Image2D;

use PDL::NiceSlice;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<3) {
    print "USE: rss2grid_new.pl INPUT_RSS position_table.txt DPIX OUTPUT_cube [OFFSET_X Y]\n";
    print "This program regularize a non regular IFU data into a datacube without";
    print "any interpolation\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$dpix=$ARGV[2];
$output_cube=$ARGV[3];
$offx=0;
$offy=0;
$no_val=1000;
if ($#ARGV==5) {
    $offx=$ARGV[4];
    $offy=$ARGV[5];
}
$offx=$offx/$dpix;
$offy=$offy/$dpix;

$n=0;
open(PT,"<$pos_table");
$line=<PT>;
chop($line);
@data=split(" ",$line);
$dpix_init=$data[1]/$dpix;
$nx_min=1e12;
$ny_min=1e12;
$nx_max=-1e12;
$ny_max=-1e12;


while($line=<PT>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $xx[$n]=$data[1]/$dpix;
    $yy[$n]=$data[2]/$dpix;
    if ($xx[$n]>$nx_max) {
	$nx_max=$xx[$n];
    }
    if ($yy[$n]>$ny_max) {
	$ny_max=$yy[$n];
    }
    if ($xx[$n]<$nx_min) {
	$nx_min=$xx[$n];
    }
    if ($yy[$n]<$ny_min) {
	$ny_min=$yy[$n];
    }
    
    $n++;
}
close(PT);

$nx=intval($nx_max-$nx_min+1);
$ny=intval($ny_max-$ny_min+1);
for ($j=0;$j<$n;$j++) {
    $xx[$j]=$xx[$j]+$offx-$nx_min;
    $yy[$j]=$yy[$j]+$offy-$ny_min;
    $ii[$j]=intval($xx[$j]);
    $jj[$j]=intval($yy[$j]);
}
print "We will create a cube of $nx X $ny pixels\n";

$input=rfits($input_rss);
($nx_rss,$ny_rss)=$input->dims;
if ($n!=$ny_rss) {
    print "The number of entries in the position table ($n)\n";
    print "does not correspond with the Y-axis of the RSS ($ny_rss)\n";
    exit;
}

$cube=zeroes($nx,$ny,$nx_rss);
$nadd=zeroes($nx,$ny);
for ($j=0;$j<$ny_rss;$j++) {
    $ii_pos_min=$ii[$j]-intval($dpix_init*1.5);
    $jj_pos_min=$jj[$j]-intval($dpix_init*1.5);
    $ii_pos_max=$ii[$j]+intval($dpix_init*1.5);
    $jj_pos_max=$jj[$j]+intval($dpix_init*1.5);
    for ($jj_pos=$jj_pos_min;$jj_pos<$jj_pos_max;$jj_pos++) {
	for ($ii_pos=$ii_pos_min;$ii_pos<$ii_pos_max;$ii_pos++) {
#	    $dist=sqrt(($ii_pos-$ii[$j])**2+($jj_pos-$jj[$j])**2);
	    $dist=sqrt(($ii_pos-$xx[$j])**2+($jj_pos-$yy[$j])**2);
#	    if ($dist<=($dpix_init*1.1)) {
	    if ($dist<=$dpix_init) {
		if (($ii_pos>0)&&($ii_pos<$nx)&&($jj_pos>0)&&($jj_pos<$ny)) {
		    my $data=$input->slice(":,$j");
		    for ($k=0;$k<$nx_rss;$k++) {
			$t=$cube->at($ii_pos,$jj_pos,$k);
			$t=$t+$input->at($k,$j);
			set($cube,($ii_pos,$jj_pos,$k),$t);
		    }
		    my $add=$nadd->slice("$ii_pos,$jj_pos"); $add .=$add+1;
		}
	    } else {
		if (($dist>$dpix_init)&&($dist<($dpix_init+$dpix))) {
		    $w=($dist-$dpix_init+0.5*$dpix)/$dpix;
		    if (($ii_pos>0)&&($ii_pos<$nx)&&($jj_pos>0)&&($jj_pos<$ny)) {
			my $data=$input->slice(":,$j");
			for ($k=0;$k<$nx_rss;$k++) {
			    $t=$cube->at($ii_pos,$jj_pos,$k);			    
			    $t=$t+($input->at($k,$j))*$w;
			    set($cube,($ii_pos,$jj_pos,$k),$t);
			}
			my $add=$nadd->slice("$ii_pos,$jj_pos"); $add .=$add+$w;
		    }
		}
	    }
	}
    }
    print "$j/$ny_rss\n";
}

for ($j=0;$j<8;$j++) {
    $F[$j]=0;
}

for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$t=$nadd->at($i,$j);
	if ($t>0) {
	    my $val=$cube->slice("$i,$j,:"); $val .=$val/$t;
	} else {
	    my $val=$cube->slice("$i,$j,:"); $val .=-$no_val;
	}
	if ($t==0) {
	    $F[0]++;
	}
	if ($t<1) {
	    $F[1]++;
	}
	if ($t==1) {
	    $F[2]++;
	}
	if (($t>1)&&($t<2)) {
	    $F[3]++;
	}
	if ($t==2) {
	    $F[4]++;
	}
	if (($t>2)&&($t<3)) {
	    $F[5]++;
	}
	if ($t==3) {
	    $F[6]++;
	}
	if ($t>3) {
	    print "$t\n";
	    $F[7]++;
	}
    }
}
for ($j=0;$j<8;$j++) {
    $F[$j]=int(($F[$j]/($nx*$ny))*10000)/100;
}
print "-------------\n";
print "N=0,   $F[0]\n";
print "0<N<1, $F[1]\n";
print "N=1,   $F[2]\n";
print "1<N<2, $F[3]\n";
print "N=2,   $F[4]\n";
print "2<N<3, $F[5]\n";
print "N=3,   $F[6]\n";
print "N>3,    $F[7]\n";
print "-------------\n";

$start_w=$input->hdr->{CRVAL1};
$delta_w=$input->hdr->{CDELT1};
$crpix1=$input->hdr->{CRPIX1};
$cube->hdr->{CRPIX1}=1;
$cube->hdr->{CRVAL1}=($nx_min+4*$dpix_init)*$dpix;
$cube->hdr->{CDELT1}=$dpix;
$cube->hdr->{CRPIX2}=1;
$cube->hdr->{CRVAL2}=($ny_max+4*$dpix_init)*$dpix;
$cube->hdr->{CDELT2}=$dpix;
$cube->hdr->{CRPIX3}=1;
$cube->hdr->{CRVAL3}=$start_w-($crpix1-1)*$delta_w;
$cube->hdr->{CDELT3}=$delta_w;

$cube->wfits($output_cube);

$nadd->wfits("nadd_new.fits");

exit;

sub intval {
    my $a=@_[0];
    my $ia=int($a);
    my $d=$a-$ia;
    if ($d>0.5) {
	$ia=$ia+1;
    }
    return $ia;
}
