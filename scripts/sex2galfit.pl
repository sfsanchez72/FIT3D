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
@param=("name","x","y","mag","re","nsersic","ab","pa");
if ($#ARGV<3) {
    print "USE: sex2galfit.pl SExtractor.cat galfit_ref fix_param_file galfit_out\n";
    print "This program transform a SEXTRACTED catalogue to a GALFIT input\n";
    print "FIX_PARAM_FILE: X Y RE NSERSIC A/B PA SUBSTRACT(0=YES)\n";
    print "SEX.CAT:\n";
    print "#   1 NUMBER          Running object number\n";
    print "#   2 X_IMAGE         Object position along x                         [pixel]\n";
    print "#   3 Y_IMAGE         Object position along y                         [pixel]\n";
    print "#   4 ALPHA_J2000     Right ascension of barycenter (J2000)           [deg]\n";
    print "#   5 DELTA_J2000     Declination of barycenter (J2000)               [deg]\n";
    print "#   6 THETA_IMAGE     Position angle (CCW/x)                          [deg]\n";
    print "#   7 ELLIPTICITY     1 - B_IMAGE/A_IMAGE\n";
    print "#   8 A_IMAGE         Profile RMS along major axis                    [pixel]\n";
    print "#   9 FWHM_IMAGE      FWHM assuming a gaussian core                   [pixel]\n";
    print "#  10 CLASS_STAR      S/G classifier output\n";
    print "#  11 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR                 [mag]\n";
    print "#  12 MAGERR_BEST     RMS error for MAG_BEST                          [mag]\n";
    print "#  13 FLAGS           Extraction flags\n";
    exit;
}

$sex_cat=$ARGV[0];
$galfit_conf=$ARGV[1];
$fix_param_file=$ARGV[2];
$galfit_out=$ARGV[3];

#
# We create a tmp galfit:
#
print "Reading $galfit_conf\n";
$model_number=0;
$nh=0;
open(FH,"<$galfit_conf");
$line="INPUT";
while ($line !~ /S\)/) {
    $line=<FH>;
    chop($line);
    $header[$nh]=$line;
    $nh++;
}
close(FH);
open(FH,"<$fix_param_file");
$line=<FH>;
chop($line);
@a_fit=split(" ",$line);
close(FH);
print "@a_fit\n";

$n=0;
open(FH,"<$sex_cat");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$id[$n]=$data[0];
	$x[$n]=$data[1];
	$y[$n]=$data[2];
	$theta[$n]=$data[5]-90;
	$ellip[$n]=$data[6];
	$ab[$n]=1-$ellip[$n];
	$r50[$n]=$data[7];
	$fwhm[$n]=$data[8];
	$star[$n]=$data[9];
	$mag[$n]=$data[16];
	$e_mag[$n]=$data[16];
	$n++;
    }
}
close(FH);
print "$sex_cat has $n entries\n";

open(FH,">$galfit_out");
for ($i=0;$i<$nh;$i++) {
    print FH "$header[$i]\n";
}
for ($i=0;$i<$n;$i++) {
    $ii=$i+1;
    print FH "\n";
    print FH "#Object number: $ii\n";
    print FH " 0) sersic # Object type\n";
	print FH " 1) $x[$i] $y[$i] $a_fit[1] $a_fit[2] # ($param[1],$param[2])\n";
	print FH " 3) $mag[$i] $a_fit[3] # ($param[3])\n";
	print FH " 4) $r50[$i] $a_fit[4] # ($param[4])\n";
	print FH " 5) 2.5 $a_fit[5] # ($param[5])\n";
	print FH " 6) 0.000      0       #     ----- \n";
	print FH " 7) 0.000      0       #     ----- \n";
	print FH " 8) $ab[$i] $a_fit[6] # ($param[6])\n";
	print FH " 9) $theta[$i] $a_fit[7] # ($param[7])\n";
	print FH "10) 0.000      0       #  diskiness(-)/boxiness(+)\n";
	print FH " Z) $a_fit[8]        #  Output option (0 = residual, 1 = Don't subtract) \n";
    }
close(FH);
exit;
