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
    print "USE: rss2cube.pl input_RSS.fits output_cube.fits nx ny\n";
    exit;
}

$input_rss=$ARGV[0];
$output_cube=$ARGV[1];
$nx=$ARGV[2];
$ny=$ARGV[3];


$nax=read_naxes($input_rss);   
@naxis=@$nax;
@a_rss=read_img($input_rss);
$check=$nx*$ny;
($start_w,$delta_w)=read_img_headers($input_rss,["CRVAL1","CDELT1"]);
if ($check!=$naxis[1]) {
    print "Error: The Y-axis of the RSS has $naxis[1] values\n";
    print "It does not correspond with $nx X $ny\n";
    exit;
}

for ($k=0;$k<$naxis[0];$k++) {    
    for ($j=0;$j<$ny;$j++) {
	for ($i=0;$i<$nx;$i++) {
#	    $a_cube[$k][$j][$i]=$a_rss[$j+($nx-1-$i)*$ny][$k];	    
	    $a_cube[$k][$i][$j]=$a_rss[$j+($nx-1-$i)*$ny][$k];	    
	}
    }
}
#$nax=[$nx,$ny,$naxis[0]];
$nax=[$ny,$nx,$naxis[0]];
system("rm $output_cube");
write_fits($output_cube,$nax,3,\@a_cube);
write_crval_cube($output_cube,[1,0,1,1,0,1,1,$start_w,$delta_w]);
exit;
