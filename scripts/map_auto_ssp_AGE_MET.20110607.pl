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
use PDL::Image2D;
use  PDL::Fit::Linfit;



use PDL::Core;
use PDL::Basic;
use PDL::Exporter;
@ISA    = qw( PDL::Exporter );
use PDL::Options ':Func';
use PDL::Slatec; # For matinv()


$vel_light=299792.458;


if ($#ARGV<1) {
    print "USE: map_auto_ssp.pl AUTO.OUT PREFIX_OUT\n";
    exit;
}

$infile=$ARGV[0];
$prefix=$ARGV[1];

$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $nx=$data[0];
    $ny=$data[1];
    $ix=$data[2];
    $iy=$data[3];
#    ($nx,$ny,$ix,$iy)=split(" ",$line);
#    $line=<FH>;
#    chop($line);
#    $nmod=$line;
    if ($n==0) {
	$pdl_chi=zeroes($nx,$ny);
	$pdl_age=zeroes($nx,$ny);
	$pdl_met=zeroes($nx,$ny);
	$pdl_Av=zeroes($nx,$ny);
	$pdl_disp=zeroes($nx,$ny);
	$pdl_vel=zeroes($nx,$ny);
	$pdl_flux=zeroes($nx,$ny);
    }
    set($pdl_chi,$ix,$iy,$data[4]);
    set($pdl_age,$ix,$iy,$data[5]);
    set($pdl_met,$ix,$iy,$data[6]);
    set($pdl_Av,$ix,$iy,$data[7]);
    set($pdl_vel,$ix,$iy,$data[8]);
    set($pdl_disp,$ix,$iy,$data[9]);
    set($pdl_flux,$ix,$iy,$data[12]);

    $n++;
}
close(FH);
$chi_name=$prefix."_chi.fits";
$pdl_chi->wfits($chi_name);

$age_name=$prefix."_age.fits";
$pdl_age->wfits($age_name);

$met_name=$prefix."_met.fits";
$pdl_met->wfits($met_name);

$Av_name=$prefix."_Av_ssp.fits";
$pdl_Av->wfits($Av_name);

$disp_name=$prefix."_disp_ssp.fits";
$pdl_disp->wfits($disp_name);

$pdl_vel=$pdl_vel*300000;
$vel_name=$prefix."_vel_ssp.fits";
$pdl_vel->wfits($vel_name);

$flux_name=$prefix."_flux_ssp.fits";
$pdl_flux->wfits($flux_name);


exit;





