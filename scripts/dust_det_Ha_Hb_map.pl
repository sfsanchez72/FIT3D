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

if ($#ARGV<2) {
    print "USE: dust_det_Ha_Hb_map.pl ratio_obs_Ha_Hb.fits ratio_expected_Ha_Hb[VALUE] AV_map.fits\n";
    exit;
}

$Rv=5.5;
$rat_obs_file=$ARGV[0];
$rat_exp=$ARGV[1];
$out=$ARGV[2];
$lambda1=4861.32;
$lambda2=6562.68;

$rat_obs=rfits($rat_obs_file);
#$rat_exp=rfits($rat_exp_file);
($nx,$ny)=$rat_obs->dims;
$Av=zeroes($nx,$ny);

$a_1=A_l(5.5,$lambda1);
$a_2=A_l(5.5,$lambda2);

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
