#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#


use Statistics::OLS;
use Math::Stat;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";




if ($#ARGV<2) {
    print "USE: gridslice.pl slice.txt fitsfile outslice\n";
    exit;
}

$slice=$ARGV[0];
$infile=$ARGV[1];
$outfile=$ARGV[2];


$ns=0;
$x_min=10e10;
$x_max=-10e10;
$y_min=10e10;
$y_max=-10e10;
$dpix=10e10;
open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	if ($x_min>$x[$ns]) {
	    $x_min=$x[$ns];
	}
	if ($y_min>$y[$ns]) {
	    $y_min=$y[$ns];
	}
	if ($x_max<$x[$ns]) {
	    $x_max=$x[$ns];
	}
	if ($y_max<$y[$ns]) {
	    $y_max=$y[$ns];
	}
	if ($ns>0) {
	    if (($dpix>abs($x[$ns]-$x[$ns-1]))&&(abs($x[$ns]-$x[$ns-1])!=0)) {
		$dpix=abs($x[$ns]-$x[$ns-1]);
	    }
	}
	$ns++;
    }
}
close(FH);

$naxes=read_naxes($infile);
($nx,$ny)=@$naxes;
@af=read_img($infile);
print "GRID=$x_min $x_max $y_min $y_max $dpix [$nx,$ny]\n";

$kk=0;
open(FH,">$outfile");
for ($k=0;$k<$ns;$k++) {
    $i=abs($x[$k]-$x_min)/$dpix;
    $j=abs($y[$k]-$y_min)/$dpix;
    $flux[$k]=$af[$j][$i];
    $kk=$k+1;
    print FH "$kk $x[$k] $y[$k] $flux[$k]\n";
}
close(FH);



#system("rm $outfile");
#write_img($outfile,$nx,$ny,\@af);
#write_wcs($outfile,1,$x_min,'ARCSEC',$dpix,0,1,$y_min,'ARCSEC',$dpix,0);

exit;


sub int_poly {
    my $point=$_[0];
    my $wave=$_[1];
    my $suma=0;
    my @data=@$point;
    my $j,$k;
    for ($j=1;$j<$#data;$j=$j+2) {
	$k=($j-1)/2;
	$suma+=$data[$j]*($wave**$k);
#	print "WAVE=$wave $suma $data[$j] $k\n";
    }
    return $suma;
}
