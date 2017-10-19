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




if ($#ARGV<4) {
    print "USE: slicegrid_grid.pl slice.txt fitsfile nx ny dpix\n";
    exit;
}

$slice=$ARGV[0];
$outfile=$ARGV[1];
$nx=$ARGV[2];
$ny=$ARGV[3];
$dpix=$ARGV[4];

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
	$flux[$ns]=$data[3];
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

#$nx=abs($x_max-$x_min)/$dpix+1;
#$ny=abs($y_max-$y_min)/$dpix+1;
print "GRID=$x_min $x_max $y_min $y_max $dpix [$nx,$ny]\n";
#exit;
$pdl=zeroes($nx,$ny);
for ($k=0;$k<$ns;$k++) {
    $i=abs($x[$k]-$x_min)/$dpix;
    $j=abs($y[$k]-$y_min)/$dpix;
    $af[$j][$i]=$flux[$k];
    set($pdl,$i,$j,$af[$j][$i]);
#   print "$i $j $af[$i][$j]\n";
}

$pdl->wfits($outfile);


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
