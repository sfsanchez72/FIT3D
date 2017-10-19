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
    print "USE: cube_shift_delta.pl input_cube.fits output_cube.fits delta_x delta_y\n";
    exit;
}

$input_cube=$ARGV[0];
$output_cube=$ARGV[1];
$dx=$ARGV[2];
$dy=$ARGV[3];


$na=read_naxes($input_cube);   
@naxis=@$na;


system("cp $input_cube $output_cube");
print "Shifting $input_cube to $output_cube\n";

for ($k=0;$k<$naxis[2];$k++) {
    $kk=$k+1;
    system("rm cube_shift_delta.fits");
    open (TMP,">cube_shift_delta.cl");
    print TMP "\n";
    $line="imshift ".$input_cube."[*,*,".$kk."] cube_shift_delta.fits ".$dx." ".$dy;
    print TMP "$line\n";
    $line="imcopy cube_shift_delta.fits ".$output_cube."[*,*,".$kk."]";
    print TMP "$line\n";
    print TMP "logout\n";
    close(TMP);
    system("/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh < cube_shift_delta.cl > cube_shift.log");
    if ($k==50*int($k/50)) {
	print "$k\n";
    }
#    print "$k, ";
}
print "DONE\n";

exit;
