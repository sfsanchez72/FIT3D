#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
#

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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<1) {
    print "USE: spec2D_plot.pl INPUT_FILE.FITS NY DEVICE [MIN MAX] [WMIN WMAX]\n [REF_LINES CUT]";
    exit;
}

$input=$ARGV[0];
$NY0=$ARGV[1];
$dev=$ARGV[2];
$def=0;
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}





$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;
$crval=$pdl->hdr->{"VAL1_1"};
$cdelt=$pdl->hdr->{"DLT1_1"};
$crpix=$pdl->hdr->{"PIX1_1"};
if ($cdelt==0) {
    $cdelt=1;
}

for ($i=0;$i<$nx;$i++) {
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
    $NY=int($ny/2);
    $flux[$i]=$pdl->at($i,$NY);
    if ($flux[$i]>$y_max) {
	$y_max=$flux[$i];
    }
    if ($flux[$i]<$y_min) {
	$y_min=$flux[$i];
    }
#    print "$wave[$i] $flux[$i]\n";
}


$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}

if ($#ARGV==6) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $def=1;
}
$ref_def=0;
if ($#ARGV==8) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $ref_lines=$ARGV[7];
    $cut_lines=$ARGV[8];
    $ref_def=1;
    $def=1;
}

if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}

#print "$y_min $y_max\n";

if ($ref_def==1) {
    $nl=0;
    open(FH,"<$ref_lines");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($data[0]>$cut_lines) {
	    $wave_line[$nl]=$data[1];
	    if ($nl>0) {
		if (abs($wave_line[$nl]-$wave_line[$nl-1])>5) {
		    $nl++;
		    }
	    } else {
		$nl++;
	    }
	}
    }
    close(FH);
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength","Flux","");
$color=1;
for ($NY=0;$NY<$ny;$NY++) {
  #    print "$NY/$ny\n";
      $k=$NY+1;    
      $val="VAL1_".$k;
     $crval=$pdl->hdr->{$val};
     $val="DLT1_".$k;
     $cdelt=$pdl->hdr->{$val};    
    my @wave;
    my @flux;
    for ($i=0;$i<$nx;$i++) {
	$wave[$i]=$crval+$cdelt*$i;
	$flux[$i]=$pdl->at($i,$NY);
	if ($NY0==-1) {
	    $n1=int($nx*0.45);
	    $n2=int($nx*0.55);
	    my $sec=$pdl->slice("$n1:$n2,$NY");
	    $norm=sumover($sec)/($n2-$n1);
	    $flux[$i]=$flux[$i]/$norm->at(0);
	}
	if ($NY0==-2) {
	    $n1=int($nx*0.45);
	    $n2=int($nx*0.55);
	    my $sec=$pdl->slice("$n1:$n2,$NY");
	    $norm=sumover($sec)/($n2-$n1);
	    $flux[$i]=$flux[$i]/$norm->at(0);
	    $flux[$i]=$flux[$i]+($y_max-$y_min)*0.05*$NY;
	}
#	print "$i $wave[$i] $flux[$i]\n";
    }
    pgsci($color);
    pgline($nx,\@wave,\@flux);
    $color++;
    if ($color>15) {
	$color=1;
    }
}
pgsch(0.9);
pgsci(2);
for ($j=0;$j<$nl;$j++) {
    if (($wave_line[$j]>$wmin)&&($wave_line[$j]<$wmax)) {
	if ($j==(2*int($j/2))) {
	    pgptxt($wave_line[$j],0.8*$y_max,90,0,"$wave_line[$j]");
	} else {
	    pgptxt($wave_line[$j],0.7*$y_max,90,0,"$wave_line[$j]");
	}
    }
}
pgclose;
pgend;

exit;
