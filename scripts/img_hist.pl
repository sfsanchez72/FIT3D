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
    print "USE: table_stats.pl FITSFILE.FITS MIN DEV NBIN\n";
    exit;
}

$infile=$ARGV[0];
$min=$ARGV[1];
$dev=$ARGV[2];
$nbin=$ARGV[3];
$pdl=rfits($infile);
($nx,$ny)=$pdl->dims;
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val_now=$pdl->at($i,$j);
	if ($val_now>$min) {
	    $val[$n]=$val_now;
	    $n++;
	}
    }
}

$slice=pdl(@val);
($mean,$rms,$median,$min,$max) = stats($slice);

print "$mean $rms $median $min $max\n";

$N=$n;
@a=list($slice);
$NN=$n/int($nbin/6+1);


pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgenv($median-3*$rms,$median+3*$rms,0,$NN,0,0);
pglabel("$infile Value","N values","$median+-$rms");
pgsfs(3);
pghist($N,\@a,$median-3*$rms,$median+3*$rms,$nbin,3);
pgsfs(2);
pghist($N,\@a,$median-3*$rms,$median+3*$rms,$nbin,3);
pgclose;
pgend;
exit;





