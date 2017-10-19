#!/usr/bin/perl
#
#
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

if ($#ARGV<5) {
    print "USE: dust_det_Ha_Hb_map.pl RV(3.1) ratio_obs_2_1.fits ratio_expected_3_1[VALUE] lambda1 lambda2 AV_map.fits\n";
    exit;
}

$Rv=$ARGV[0];
$rat_obs_file=$ARGV[1];
$rat_exp=$ARGV[2];
$lambda1=$ARGV[3];
$lambda2=$ARGV[4];
$out=$ARGV[5];


$rat_obs=rfits($rat_obs_file);
#$rat_exp=rfits($rat_exp_file);
($nx,$ny)=$rat_obs->dims;
$Av=zeroes($nx,$ny);

$a_1=A_l($Rv,$lambda1);
$a_2=A_l($Rv,$lambda2);

for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$rat_obs_val=$rat_obs->at($i,$j);
	$rat_exp_val=$rat_exp;#->at($i,$j);
	if ($rat_obs_val>$rat_exp_val) {
	    $Av_val=(2.5*log10($rat_obs_val/$rat_exp_val))/($a_1-$a_2);
	    set($Av,$i,$j,$Av_val);
	}
    }
}
#print "$Av\n";
$Av->wfits($out);

exit;
