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

@x=(10,20,30,40,50,100,200,300,400,500,1000,2000,3000,4000,5000,10000,20000,30000,100000);
@y=(1.43,1.41,1.40,1.39,1.37,1.33,1.25,1.15,1.10,1.06,0.85,0.73,0.60,0.58,0.53,0.49,0.47,0.46,0.45);

for ($i=0;$i<($#x+1);$i++) {
    $lx[$i]=log10($x[$i]);
}

if ($#ARGV<2) {
    print "USE: SII_ratio.pl SII_ratio Temperature(1e4) DEV\n";
    exit;
}
$ratio=$ARGV[0];
$T=$ARGV[1];
$dev=$ARGV[2];


$pdl_lx=pdl(@lx);
$pdl_y=pdl(@y);
$npoly=8;
($pdl_yp,$coeff) = fitpoly1d $pdl_lx,$pdl_y,$npoly;
#$npoly=8;
#($pdl_xp,$coeff2) = fitpoly1d $pdl_yp,$pdl_lx,$npoly;
for ($j=0;$j<$npoly;$j++) {
    $c[$j]=$coeff->at($j);
}
@yp=list($pdl_yp);


$min=1e12;
for ($i=100;$i<500;$i++) {
    $lden=$i/100;
    $rat=0;
    for ($j=0;$j<$npoly;$j++) {
	$rat=$rat+$c[$j]*($lden**($j));
    }
    $dist=abs($rat-$ratio);
    if ($dist<$min) {
	$min=$dist;
	$LDENSITY=$lden;
    }
}





#@xp=list($pdl_xp);
#for ($i=0;$i<($#x+1);$i++) {
#    print "$i $lx[$i] $y[$i] $yp[$i]\n";
#}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv(0.99,5.01,0.0,1.6,10,0);
pglabel("Ne(10\\u4\\d/T)\\u1/2","[SII]6716/6731","");
pgpoint(($#x+1),\@lx,\@y,3);
pgsci(2);
pgline(($#x+1),\@lx,\@yp);
pgsci(3);
pgpoint(1,[$LDENSITY],[$ratio],5);
#pgline(($#x+1),\@xp,\@y);
pgclose;
pgend;

$DENSITY=10**($LDENSITY);
$FACT=(10**2)/sqrt($T);
$REAL_DENSITY=$DENSITY/$FACT;
print "$REAL_DENSITY\n";
#print "$FACT";

exit;
