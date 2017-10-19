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


if ($#ARGV<2) {
    print "USE: table_hist_map.pl CUBE.fits DEV NBIN [GREP]\n";
    exit;
}

$infile=$ARGV[0];
$dev=$ARGV[1];
$nbin=$ARGV[2];
$no_grep=0;
if ($#ARGV==3) {
    $no_grep=1;
    $grep=$ARGV[3];
}
$n=0;
$pdl=rfits($infile);
($nx,$ny)=$pdl->dims;

$n=0;
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$data=$pdl->at($i,$j);
	if ($data>0) {
	    $val[$n]=$data;
	    $n++;
	}
    }
}

$slice=pdl(@val);

#$slice=rfits($infile);
($mean,$rms,$median,$min,$max) = stats($slice);

#@val=list($slice);


print "$mean $rms $median $min $max\n";

$N=$n;
@a=list($slice);
$NN=$n/int($nbin/30+1);

$median=median(@a);
$rms=sigma(@a);

$k=0;
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$data=$pdl->at($i,$j);
	if ($data<($median+$rms)) {
	    $valk[$k]=$data;
	    $k++;
	}
    }
}

$mediank=median(@valk);
$rmsk=sigma(@valk);
$cut=$mediank+3*$rmsk;
print "$mediank $rmsk $cut\n";
pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($median-5*$rms,$median+5*$rms,0,$NN,0,0);
pglabel("$infile Value","N values","$median+-$rms");
pgsfs(3);
pghist($N,\@a,$median-5*$rms,$median+5*$rms,$nbin,3);
pgsfs(2);
pghist($N,\@a,$median-5*$rms,$median+5*$rms,$nbin,3);
pgsci(2);
pgline(2,[$median+$rms,$median+$rms],[0,$n]);
pgsci(3);
pgline(2,[$mediank+3*$rmsk,$mediank+3*$rmsk],[0,$n]);
pgclose;
pgend;
exit;





