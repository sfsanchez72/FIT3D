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


if ($#ARGV<3) {
    print "USE: correct_galExt_rss.pl INPUT_FILE.FITS A_lambda wavelength OUTPUT_FILE.FITS\n";
    exit;
}

$input=$ARGV[0];
$Al=$ARGV[1];
$lambda=$ARGV[2];
$output=$ARGV[3];
$def=0;

$pdl=rfits($input);
$h=$pdl->gethdr;
($nx,$ny)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL1"};
$cdelt=$pdl->hdr->{"CDELT1"};
$crpix=$pdl->hdr->{"CRPIX1"};

$a_V=A_l(3.1,5505);
$a_l=A_l(3.1,$lambda);

$AV=($a_V/$a_l)*$Al;

$pdl_C=zeroes($nx);

for ($i=0;$i<$nx;$i++) {
    $wave=$crval+$cdelt*$i;
    $C=$AV*A_l(3.1,$wave);
    $F=(10**(0.4*($C)));
    set($pdl_C,$i,$F);
}
for ($j=0;$j<$ny;$j++) {
    my $flux=$pdl->slice(":,($j)");
    $flux .= $flux*$pdl_C;
}

$pdl->wfits($output);



exit;
