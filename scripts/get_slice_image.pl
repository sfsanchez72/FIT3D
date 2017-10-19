#!/usr/bin/perl
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

if ($#ARGV<6) {
    print "USE: get_slice_image.pl INPUT_IMAGE.fits POS_TABLE.txt PIX_SCALE DELTA_X_PIX DELTA_Y_PIX output_slice.txt FACTOR [DX DY]\n";
    exit;
}

$input_image=$ARGV[0];
$slice=$ARGV[1];
$pix_scale=$ARGV[2];
$delta_x=$ARGV[3];
$delta_y=$ARGV[4];
$out_file=$ARGV[5];
$factor=$ARGV[6];
$DX=$ARGV[7];
$DY=$ARGV[8];


$ns=0;
open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	$ns++;
    } else {
	$header=$line;
	$type=$data[0];
	$aper=$data[1];
    }
}
close(FH);



$b=rfits("$input_image");
$h=$b->gethdr;
($nx,$ny)=$b->dims;
$PHOTFLAM=$b->hdr->{PHOTFLAM};
$ZPOINT=$b->hdr->{ZPOINT};


$rat=$factor*$PHOTFLAM/$ZPOINT;




for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {

    }
}
print "$nx,$ny $delta_x,$delta_y";
open(OUT,">$out_file");
print OUT "$header\n";
$sum=0;
$mag=0;
for ($j=0;$j<$ns;$j++) {
    $xc=$x[$j]-2+$DX;
    $yc=$y[$j]-2+$DY;
    $I=-$yc/$pix_scale+($ny-$delta_x)-26/$pix_scale;
    $J=($xc/$pix_scale)-($nx-$delta_y);#-$xc/$pix_scale;#-$delta_x;
    $val=0;
    if ($type eq "C") {
	$i_min=int($I-4*$aper);
	$i_max=int($I+4*$aper);
	$j_min=int($J-4*$aper);
	$j_max=int($J+4*$aper);
	for ($ii=$i_min;$ii<$i_max;$ii++) {
	    for ($jj=$j_min;$jj<$j_max;$jj++) {
		$flux=$b->at($ii,$jj);
		$rad=sqrt(($ii-$I)**2+($jj-$J)**2);
		$rad=$rad*$pix_scale;
		if ($rad<=$aper) {
		    $val=$val+$flux;
		}
	    }
	}
    }
    if ($type eq "S") {
	$i_min=int($I-$aper);
	$i_max=int($I+$aper);
	$j_min=int($J-$aper);
	$j_max=int($J+$aper);
	for ($ii=$i_min;$ii<$i_max;$ii++) {
	    for ($jj=$j_min;$jj<$j_max;$jj++) {
		$flux=$b->at($ii,$jj);
		$val=$val+$flux;
	    }
	}
    }
    $val=$val*$rat;
    $sum=$sum+$val;
    print OUT "$id[$j] $x[$j] $y[$j] $val\n";
}
close(OUT);

$mag=2.5*log10($rat/$sum*1e16);
print "MAG = $mag\n";


exit;






