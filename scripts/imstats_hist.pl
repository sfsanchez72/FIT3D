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


if ($#ARGV<6) {
    print "USE: imstats_hist.pl INPUT.FITS X_C Y_C X_WIDTH Y_WIDTH DEV NBIN\n";
    exit;
}

$infile=$ARGV[0];
$x_c=$ARGV[1];
$y_c=$ARGV[2];
$dx=$ARGV[3];
$dy=$ARGV[4];
$dev=$ARGV[5];
$nbin=$ARGV[6];
$x_ini=$x_c-$dx;
$x_end=$x_c+$dx;
$y_ini=$y_c-$dy;
$y_end=$y_c+$dy;
$pdl=rfits($infile);
($nx,$ny)=$pdl->dims;
$slice=$pdl->slice("$x_ini:$x_end,$y_ini:$y_end");
($mean,$rms,$median,$min,$max) = stats($slice);

print "$mean $rms $median $min $max\n";

$N=($x_end-$x_ini+1)*($y_end-$y_ini+1);
@a=list($slice);
$NN=$nx*$ny/7;


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





