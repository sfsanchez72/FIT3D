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

if ($#ARGV<5) {
    print "USE: auto_ssp_elines_several_Av_log_rss_loop.pl SPEC1.RSS.fits LIST_OF_SSP.txt PREFIX.OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
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
$call="date > date.start";
system($call);


$unc_file=$ARGV[0];
$clean_file="clean_".$ARGV[0];


$back_list=$ARGV[1];
$outfile=$ARGV[2];
$out_elines="elines_".$outfile;
$out_single="single_".$outfile;
$out_fit="fit_".$outfile;
$out_coeffs="coeffs_".$outfile;




$mask_list=$ARGV[3];
$config_file=$ARGV[4];
$plot=$ARGV[5];

if ($plot==2) {
    $plot=1;
    $dev_plot=$outfile.".ps/CPS";
    $dev_plot_single="single_".$outfile.".ps/CPS";
} else {
    $dev_plot="/xs";
    $dev_plot_single="/xs";
}

$smooth=1;
$MIN_CHISQ=1e12;


$out_file="junk.junk";
$factor=1;
$box=1;


$def=0;
if ($#ARGV==7) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $def=1;
}

if ($#ARGV==9) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $def=2;
}

if ($#ARGV==10) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $def=2;
}

$input_redshift=0;
if ($#ARGV==11) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $def=2;
}



if ($#ARGV==14) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $def=2;
}
#print "$#ARGV $d_redshift\n"; 
open(FH,"<$config_file");
$line=<FH>;
chop($line);
($redshift,$d_redshift,$min_redshift,$max_redshift)=split(" ",$line);

$line=<FH>;
chop($line);
($sigma,$d_sigma,$min_sigma,$max_sigma)=split(" ",$line);

if ($#ARGV==18) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $sigma=$ARGV[15];
    $d_sigma=$ARGV[16];
    $min_sigma=$ARGV[17];
    $max_sigma=$ARGV[18];
    $def=2;
}




$line=<FH>;
chop($line);
($Av_IN,$d_Av_IN,$min_Av,$max_Av)=split(" ",$line);
if ($#ARGV==22) {
    $min=$ARGV[6];
    $max=$ARGV[7];
    $min_wave=$ARGV[8];
    $max_wave=$ARGV[9];
    $elines_mask=$ARGV[10];
    $input_redshift=$ARGV[11];
    $input_d_redshift=$ARGV[12];
    $input_min_redshift=$ARGV[13];
    $input_max_redshift=$ARGV[14];
    $sigma=$ARGV[15];
    $d_sigma=$ARGV[16];
    $min_sigma=$ARGV[17];
    $max_sigma=$ARGV[18];
    $Av_IN=$ARGV[19];
    $d_Av_IN=$ARGV[20];
    $min_Av=$ARGV[21];
    $max_Av=$ARGV[22];
    $def=2;
}




if ($input_redshift>0) {
    $redshift=$input_redshift;
    $d_redshift=$input_d_redshift;
    $min_redshift=$input_min_redshift;
    $max_redshift=$input_max_redshift;
}
$Av_ini=$Av_IN;

if ($d_redshift!=0) {
    $fit_redshift=1;
} else {
    $fit_redshift=0;
}

print "FIT_RED $fit_redshift $d_redshift $#ARGV\n";# <stdin>;

$NJ=0;
for ($j=0;$j<=$#ARGV;$j++) {
    $entry[$j]=$ARGV[$j];
    $NJ++;
}

$pdl=rfits($unc_file);
@dims=$pdl->dims;
$NN=$dims[1];
$wmid=0.5*($min_wave+$max_wave);

for ($j=0;$j<$NN;$j++) {    
    $call="img2spec.pl ".$unc_file." ".$j." spec_loop.txt";
    system($call);
    $entry[0]="spec_loop.txt";
    $entry[2]=$outfile.".".$j.".out";
    $call="auto_ssp_elines_several_Av_log_loop.pl ";
    for ($i=0;$i<$NJ;$i++) {
	$call=$call." ".$entry[$i];
    }
    system($call);
    $call="cp mod_joint_spec.fits mod_joint_spec.".$j.".fits";
    system($call);
    $call="cp model_spec.fits model_spec.".$j.".fits";
    system($call);
    $call="cp res_joint_spec.fits res_joint_spec.".$j.".fits";
    system($call);
}
close(OUT);


$call="date > date.end";
system($call);

exit;


