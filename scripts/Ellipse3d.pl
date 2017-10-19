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

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";
$ellipse="/home/sanchez/sda2/code/tests2/Ellipse";

if ($#ARGV<9) {
    print "USE: Ellipse3d.pl INPUT_CUBE.fits MOD_CUBE.fits RES_CUBE.fits FIX X_C Y_C PA ELLIP BACK DBACK\n";
    print "GALFIT_CONF: A) e_object.fits, B) output.fits\n";
    exit;
}

$input_cube=$ARGV[0];
$mod_cube=$ARGV[1];
$res_cube=$ARGV[2];
$fix=$ARGV[3];
$x_c=$ARGV[4];
$y_c=$ARGV[5];
$pa=$ARGV[6];
$ellip=$ARGV[7];
$back=$ARGV[8];
$dback=$ARGV[9];

$pdl=rfits($input_cube,{data=>0});

$nx=$pdl->{"NAXIS1"};
$ny=$pdl->{"NAXIS2"};
$nz=$pdl->{"NAXIS3"};

#system("rm fit.log");
#
#We create the output files
#
system("rm $mod_cube");
system("rm $res_cube");
open(TMP,">tmp_e.cl");
print TMP "\n";
print TMP "\n";
print TMP "\n";
print TMP "imcopy $input_cube $mod_cube\n";
print TMP "imcopy $input_cube $res_cube\n";
print TMP "\n";
print TMP "logout\n";
print TMP "\n";
close(TMP);
$call=$cl." < tmp_e.cl > junk.tmp";
system($call);
    



#
# Loop over the Z direction:
#
for ($k=0;$k<$nz;$k++) {
    $kk=$k+1;
    print "$k/$nz...";
    system("rm e_object.fits");
    open(TMP,">tmp_e.cl");
    print TMP "\n";
    print TMP "imcopy $input_cube\[*,*,$kk\] e_object.fits\n";
    print TMP "\n";
    print TMP "logout\n";
    print TMP "\n";
    close(TMP);
    $call=$cl." < tmp_e.cl > junk.tmp";
    system($call);
    print "copied...";
    $call=$ellipse." e_object.fits ".$fix." ".$x_c." ".$y_c." ".$pa." ".$ellip." ".$back." ".$dback;
    system($call);
    print "modelled...";
    open(TMP,">tmp_e.cl");
    print TMP "\n";
    print TMP "imcopy e_mod.fits $mod_cube\[*,*,$kk\]\n";
    print TMP "\n";
    print TMP "imcopy e_res.fits $res_cube\[*,*,$kk\]\n";
    print TMP "\n";
    print TMP "logout\n";
    print TMP "\n";
#    print TMP "display output.fits\[1\] 1 zs- zr- z1=-0.001 z2=0.02\n";
#    print TMP "display output.fits\[2\] 2 zs- zr- z1=-0.001 z2=0.02\n";
#    print TMP "display output.fits\[3\] 3 zs- zr- z1=-0.001 z2=0.02\n";
    close(TMP);
    $call=$cl." < tmp_e.cl > junk.tmp";
    system($call);
    print "DONE\n";
    open(TMP,">display_e.cl");
    print "\n";
    print TMP "display e_org.fits 1 zs- zr- z1=-0.001 z2=0.02\n";
    print TMP "display e_mod.fits 2 zs- zr- z1=-0.001 z2=0.02\n";
    print TMP "display e_res.fits 3 zs- zr- z1=-0.001 z2=0.02\n";
    close(TMP);
}
#system("cp fit.log $logfile");
print "Model done\n";

exit;




