#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#f
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
$red_elines=0.0;

   

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<5) {
    print "USE: auto_ssp_elines_several_Av_log_rss.pl SPEC1.RSS.fits BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
    print "CONFIG_FILE:\n";
    print "redshift delta_redshift min_redshift max_redshift\n";
    print "sigma delta_sigma min_sigma max_sigma\n";
    print "Av delta_Av min_Av max_Av [Same range for all]\n";
    print "N_SYSTEMS\n";
    print "(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY\n";
    print "...\n";
    print "(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY\n";
    print "MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX\n";
    print "start_w_peak end_w_peak\n";
    print "wavelength_to_norm width_AA new_back_templates.fits\n";


    exit;
}

$infile=$ARGV[0];
$pdl=rfits($infile);
($nx,$ny,$nz)=$pdl->dims();

for ($j=0;$j<$#ARGV+1;$j++) {
    $PARAM[$j]=$ARGV[$j];
}
$NP=$#ARGV+1;
$PARAM[2]=" tmp.out";
$elines="elines_".$ARGV[2];
$call="rm -f ".$elines;
$coeffs="coeffs_".$ARGV[2];
$call="rm -f ".$elines;
system($call);
$call="rm -f auto.log";
system($call);
open(OUT,">$ARGV[2]");
for ($i=0;$i<$ny;$i++) {
    $call="rm -f spec.auto.txt";
    system($call);
    $call="rm -f tmp.out";
    system($call);
    $call="img2spec.pl ".$infile." ".$i." spec.auto.txt";
    system($call);
    $call="auto_ssp_elines_several_Av_log.pl spec.auto.txt ";
    for ($j=1;$j<$NP;$j++) {
	$call=$call." ".$PARAM[$j];
    }
    $call=$call." >> auto.log";
    system($call);
    open(FH,"<tmp.out");
    $line=<FH>;
    if ($i==0) {
	print OUT "$line";
    }
    $line=<FH>;
    print OUT "$line";
    close(FH);
    $call="cat elines_tmp.out >> ".$elines;
    system($call);
    $call="cat coeffs.out >> ".$coeffs;
    system($call);
    print "$i/$ny DONE\n";
}
close(OUT);

exit;
