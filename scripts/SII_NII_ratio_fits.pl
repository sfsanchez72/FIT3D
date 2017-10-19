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

if ($#ARGV<4) {
    print "USE: SII_NII_ratio.pl SII_ratio NII_ratio OUT_DENSITY OUT_TEMP DEV\n";
    exit;
}
$ratio_file=$ARGV[0];
$ratio_NII_file=$ARGV[1];
$out_den=$ARGV[2];
$out_temp=$ARGV[3];
$dev=$ARGV[4];

$ratio=rfits($ratio_file);
$ratio_NII=rfits($ratio_NII_file);



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


($nx,$ny)=$ratio->dims;

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
#pgenv(0.99,5.01,0.0,1.6,10,0);
pgenv(0.099,6.01,0.0,1.6,10,0);
pglabel("Ne(10\\u4\\d/T)\\u1/2","[SII]6716/6731","");
pgpoint(($#x+1),\@lx,\@y,3);
pgsci(2);
pgline(($#x+1),\@lx,\@yp);

$LDENSITY=zeroes($nx,$ny);
for ($ii=0;$ii<$nx;$ii++) {
    for ($jj=0;$jj<$ny;$jj++) {
	$min=1e12;
#	for ($i=100;$i<500;$i++) {
	for ($i=100;$i<500;$i++) {
	    $lden=$i/100;
	    $rat=0;
	    for ($j=0;$j<$npoly;$j++) {
		$rat=$rat+$c[$j]*($lden**($j));
	    }
	    $val_ratio=$ratio->at($ii,$jj);
	    $dist=abs($rat-$val_ratio);
	    if ($dist<$min) {
		$min=$dist;
		$val_LDENSITY=$lden;
	    }
	}

	pgsci(3);
	pgpoint(1,[$val_LDENSITY],[$val_ratio],5);


	set($LDENSITY,$ii,$jj,$val_LDENSITY);
    }
    print "$ii/$nx\n";
}

pgclose;
pgend;

$DENSITY=10**($LDENSITY);
$T=2.5e4/(log($ratio_NII*(1+(2.5e-5)*$DENSITY)/6.91));
$FACT=(10**2)/sqrt($T);
$REAL_DENSITY=$DENSITY/$FACT;

$REAL_DENSITY->wfits($out_den);
$T->wfits($out_temp);


exit;
