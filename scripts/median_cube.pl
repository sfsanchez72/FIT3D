#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#


use Statistics::OLS;
use Math::Stat;
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
    print "USE: median_cube.pl input_cube.fits output_cube.fits XBOX YBOX\n";
    exit;
}

$input_cube=$ARGV[0];
$output_cube=$ARGV[1];
$xbox=$ARGV[2];
$ybox=$ARGV[3];


$na=read_naxes($input_cube);   
@naxis=@$na;


open (TMP,">median_cube.cl");
print TMP "\n";
print TMP "!cp $input_cube $output_cube\n";
print TMP "\n";
system("rm median_cube.fits");
for ($k=0;$k<$naxis[2];$k++) {
    $kk=$k+1;
    print TMP "!rm median_cube.fits\n";
    $line="median ".$input_cube."[*,*,".$kk."] median_cube.fits ".$xbox." ".$ybox;
    print TMP "$line\n";
    $line="imcopy median_cube.fits ".$output_cube."[*,*,".$kk."]";
    print TMP "$line\n";
    print TMP "\n";
#    print "$k/$naxis[2]\n";
}
print TMP "logout\n";
close(TMP);
system("/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh < median_cube.cl > median_cube.log");


print "DONE\n";

exit;
