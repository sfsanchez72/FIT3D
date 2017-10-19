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


if ($#ARGV<5) {
    print "USE: cube_shift.pl input_cube.fits output_cube.fits centers_list.txt X0 Y0 shift.cl  [NSTART]\n";
    exit;
}

$input_cube=$ARGV[0];
$output_cube=$ARGV[1];
$centers_list=$ARGV[2];
$x0=$ARGV[3];
$y0=$ARGV[4];
$shift_cl=$ARGV[5];
$nstart=0;
if ($#ARGV==6) {
    $nstart=$ARGV[6];
}

$na=read_naxes($input_cube);   
@naxis=@$na;

$n=0;
open(FH,"<$centers_list");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $dx[$n]=$x0-$data[1];
    $dy[$n]=$y0-$data[2];
    $dx[$n]=substr($dx[$n],0,6);
    $dy[$n]=substr($dy[$n],0,6);
    $n++;
}
close(FH);



#system("cp $input_cube $output_cube");
#print "Shifting $input_cube to $output_cube\n";

open (TMP,">$shift_cl");
print TMP "\n";
print TMP "!cp $input_cube $output_cube\n";
print TMP "\n";
system("rm cube_shift.fits");
for ($k=0;$k<$naxis[2];$k++) {
    $kk=$k+1;
    $delta_x=$dx[$k+$nstart];
    $delta_y=$dy[$k+$nstart];
    print TMP "imdel cube_shift.fits\n";
    $line="imshift ".$input_cube."[*,*,".$kk."] cube_shift.fits ".$delta_x." ".$delta_y;
    print TMP "$line\n";
    $line="imcopy cube_shift.fits ".$output_cube."[*,*,".$kk."]";
    print TMP "$line\n";
    print TMP "\n";
    #print "$k,";
#    if ($k==50*int($k/50)) {
#	print "$dx[$k] $dy[$k] $k\n";
#    }
}
    print TMP "logout\n";
    close(TMP);
    system("/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh < $shift_cl > cube_shift.log");


print "DONE\n";

exit;
