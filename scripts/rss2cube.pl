#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#



use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<3) {
    print "USE: rss2cube.pl input_rss.fits output_cube.fits nx ny\n";
    exit;
}

$input_rss=$ARGV[0];
$output_cube=$ARGV[1];
$nx=$ARGV[2];
$ny=$ARGV[3];



$in_pdl=rfits($input_rss);
$h=$in_pdl->gethdr;
($nz,$npx)=$in_pdl->dims;
$out_pdl=zeroes($nx,$ny,$nz);
for ($k=0;$k<$nz;$k++) {    
    for ($j=0;$j<$ny;$j++) {
	for ($i=0;$i<$nx;$i++) {
	    $px=$j+($nx-1-$i)*$ny;
	    $val=$in_pdl->at($k,$px);
	    set($out_pdl,$i,$j,$k,$val);
	}
    }
}
$out_pdl->sethdr($h);
$out_pdl->wfits($output_cube);

exit;
