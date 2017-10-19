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
    print "USE: extract_aper_spec_cube.pl INPUT_CUBE X Y APERTURE OUTSPEC \n";
    print "This program extract 1D spectra from a cube at the location X,Y with an\n";
    print "Aperture of R=APERTURE\n";
    exit;
}

$input_cube=$ARGV[0];
$x=$ARGV[1];
$y=$ARGV[2];
$aperture=$ARGV[3];
$out_spec=$ARGV[4];


$pdl_input=rfits($input_cube);

($nx,$ny,$nz)=$input->dims;

if ($n!=$ny_rss) {
    print "The number of entries in the position table ($n)\n";
    print "does not correspond with the Y-axis of the RSS ($ny_rss)\n";
    exit;
}

$cube=zeroes($nx,$ny,$nx_rss);
$nadd=zeroes($nx,$ny);
for ($j=0;$j<$ny_rss;$j++) {
    $ii_pos_min=$ii[$j]-intval($dpix_init*1.3);
    $jj_pos_min=$jj[$j]-intval($dpix_init*1.3);
    $ii_pos_max=$ii[$j]+intval($dpix_init*1.3);
    $jj_pos_max=$jj[$j]+intval($dpix_init*1.3);
    for ($jj_pos=$jj_pos_min;$jj_pos<$jj_pos_max;$jj_pos++) {
	for ($ii_pos=$ii_pos_min;$ii_pos<$ii_pos_max;$ii_pos++) {
	    $dist=sqrt(($ii_pos-$ii[$j])**2+($jj_pos-$jj[$j])**2);
	    if ($dist<=($dpix_init*1.07)) {
#	    if ($dist<=($dpix_init*1)) {
		if (($ii_pos>0)&&($ii_pos<$nx)&&($jj_pos>0)&&($jj_pos<$ny)) {
		    my $data=$input->slice(":,$j");
		   # my $t=$cube->slice("$ii_pos,$jj_pos,:"); $t .=$t+$data;
		    for ($k=0;$k<$nx_rss;$k++) {
			$t=$cube->at($ii_pos,$jj_pos,$k);
			$t=$t+$input->at($k,$j);
			set($cube,($ii_pos,$jj_pos,$k),$t);
		    }
		    my $add=$nadd->slice("$ii_pos,$jj_pos"); $add .=$add+1;
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
	if (($t>0)&&($t<1)) {
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
$cube->hdr->{CRPIX1}=1;
$cube->hdr->{CRVAL1}=$nx_min;
$cube->hdr->{CDELT1}=$dpix_init;
$cube->hdr->{CRPIX2}=1;
$cube->hdr->{CRVAL2}=$ny_min;
$cube->hdr->{CDELT2}=$dpix_init;
$cube->hdr->{CRPIX3}=1;
$cube->hdr->{CRVAL3}=$start_w;
$cube->hdr->{CDELT3}=$delta_w;

$cube->wfits($output_cube);

$nadd->wfits("nadd.fits");
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
