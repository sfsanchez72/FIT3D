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


#$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<1) {
    print "USE: median_column.pl RSS.fits MEDIAN_COLUMN.fits\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];


$pdl=rfits($infile);
($nx,$ny)=$pdl->dims;
$out=zeroes($nx,$ny);

for ($j=0;$j<$ny;$j++) {
    my $sec=$pdl->slice(":,$j");
    my $median=median($sec);
    my $ones=ones($nx);
    my $t = $out->slice(":,$j");
    $t .= $median*$ones;
}

$out->wfits($out_file);
print "DONE\n";

exit;
