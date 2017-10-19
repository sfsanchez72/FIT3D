#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<3) {
    print "USE: plot_spots.pl file1.spots file2.spots DEV SCALE\n";
    exit;
}
$file1=$ARGV[0];
$file2=$ARGV[1];
$dev=$ARGV[2];
$scale=$ARGV[3];

$n=0;
open(FH,"<$file1");
$line=<FH>;
chop($line);
@data=split(" ",$line);
$az1=$data[1]-180;
$el1=$data[2];
$name1=$file1;
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $x1[$n]=$data[2];
    $y1[$n]=$data[3];
    $n++;
}
close(FH);

$n=0;
open(FH,"<$file2");
$line=<FH>;
chop($line);
@data=split(" ",$line);
$az2=$data[1]-180;
$el2=$data[2];
$name2=$file2;
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $x2[$n]=$data[2];
    $y2[$n]=$data[3];
    $n++;
}
close(FH);

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv(-200,1500,-200,2500,1,0);
pglabel("X pix","Y pix","$file1 vs. $file2 (X$scale)");
pgsci(1);
pgsch(0.7);
pgsci(8);
for ($i=0;$i<$n;$i++) {
    $dx[$i]=$x2[$i]-$x1[$i];
    $dy[$i]=$y2[$i]-$y1[$i];
#    print "$x1[$i] $x2[$i] $y1[$i] $y2[$i] $dx[$i] $dy[$i]\n";
#   pgarro($x1[$i],$y1[$i],$x1[$i]+$dx[$i]*$scale,$y1[$i]+$dy[$i]*$scale);

}

$x_med=median(@dx);
$x_sig=sigma(@dx);
$y_med=median(@dy);
$y_sig=sigma(@dy);

#print "$x_med $y_med $x_sig $y_sig\n";

$nb=0;
for ($i=0;$i<$n;$i++) {
    if (abs($dx[$i]-$x_med)>1.5*$x_sig) {
	$dx[$i]=$dx[$i-1];
	$nb++;
    }
    if (abs($dy[$i]-$y_med)>1.5*$y_sig) {
	$dy[$i]=$dy[$i-1];
	$nb++;
    }
#    $dy[$i]=$y2[$i]-$y1[$i];
#    print "$x1[$i] $x2[$i] $y1[$i] $y2[$i] $dx[$i] $dy[$i]\n";
    pgarro($x1[$i],$y1[$i],$x1[$i]+$dx[$i]*$scale,$y1[$i]+$dy[$i]*$scale);

}

pgsch(0.5);
pgsci(3);
pgpoint($n,\@x1,\@y1,16);
pgsci(2);
pgpoint($n,\@x2,\@y2,22);

$x_med=median(@dx);
$x_sig=sigma(@dx);
$y_med=median(@dy);
$y_sig=sigma(@dy);
$fb=apr($nb/(2*$n));

print "$x_med $y_med $x_sig $y_sig $fb\n";


pgsci(1);
pgpoint(1,[0],[0],16);
pgarro(0,0,$scale*$x_med,$scale*$y_med);
pgsch(1.7);
pgclose;
pgend;


exit;
