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

if ($#ARGV<18) {
    print "USE: auto_ssp_elines_several_Av_log_loop.pl SPEC1.txt LIST_OF_SSP.txt OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
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


$call="rm -f mod_joint_spec.txt.*.txt";
system($call);
$call="rm -f model_spec.txt.*.txt ";
system($call);
$call="rm -f res_joint_spec.txt.*.txt ";
system($call);


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



$wmid=0.5*($min_wave+$max_wave);
$nf=0;
open(FILE,"<$back_list");
while($file=<FILE>) {
    chop($file);
    if ($file !~ "#") {
	$files[$nf]=$file;
	$nn[$nf]=0;
	$k=0;
	open(NAME,"<$file");
	while($name=<NAME>) {
	    chop($name);
	    if ($name !~ "#") {
		$nam_file[$nf][$k]=$name;
		$k++;
	    }
	}
	close(NAME);
	$nn[$nf]=$k;
       	$nf++;
    }
}
close(FILE);

#print "nf=$nf\n";
$NC=0;
$name_now=loop("",0);

$call="rm -f ".$out_elines;
system($call);
$call="rm -f ".$out_coeffs;
system($call);
open(OUT,">$outfile");
for ($j=0;$j<$NC;$j++) {
    @list_now=split(" ",$LIST[$j]);
    open(SSP,">ssp.list");
    for ($jj=0;$jj<=$#list_now;$jj++) {
	print SSP "$list_now[$jj]\n";
    }
    close(SSP);
    $call="spectra2img_norm_resample.pl ssp.list ssp.fits ".$wmid." 3300 1";
    system($call);
    $entry[1]="ssp.fits";
    $entry[2]="ssp.out";
    

    #
    # We determine the redshift and sigma and fix it
    #
#    print "j = $j\n"; exit;
    if ($j==0) {
	open(CONF,">redshift.config");
	print CONF "$redshift $d_redshift $min_redshift  $max_redshift\n";
	print CONF "$sigma $d_sigma $min_sigma $max_sigma\n";
	print CONF "$Av_IN $d_Av_IN $min_Av $max_Av\n";
	print CONF "0\n";
	print CONF "0.0001  1  0.05\n";
	print CONF "6628 6760\n";
	close(CONF);
	$call="auto_ssp_elines_several_Av_log.pl ".$entry[0]." ssp.fits outfile.out none redshift.config ".$plot;
#	print "$call\n"; exit;
	system($call);
	open(RED,"<outfile.out");
	$lin=<RED>;
	$lin=<RED>;
	chop($lin);
	@a_lin=split(" ",$lin);
	$redshift=$a_lin[4];
	$sigma=$a_lin[5];
	$entry[11]=$redshift;
	$entry[12]=0;
	$entry[13]=0.99*$redshift;
	$entry[14]=1.01*$redshift;
	$entry[15]=$sigma;
	$entry[16]=0;
	$entry[17]=0.99*$sigma;
	$entry[18]=1.01*$sigma;
	print "REDSHIFT = $redshift SIGMA = $sigma\n"; 
#exit;
    }



    $call="auto_ssp_elines_several_Av_log.pl ";
    for ($i=0;$i<$NJ;$i++) {
	$call=$call." ".$entry[$i];
    }
#    print "$call\n";
    system($call);
#    print "Press Enter\n"; <stdin>;
    open(FH,"<ssp.out");
    $out_now=<FH>;
    if ($j==0) {
	chop($out_now);
	print OUT "$out_now\n";
    }
    $out_new=<FH>;
    close(FH);
    chop($out_new);
    print OUT "$out_new\n";
    $call="cat elines_ssp.out >> ".$out_elines;
    system($call);
    $call="cat coeffs.out >> ".$out_coeffs;
    system($call);
    $call="cp mod_joint_spec.txt mod_joint_spec.txt.".$j.".txt";
    system($call);
    $call="cp model_spec.txt model_spec.txt.".$j.".txt";
    system($call);
    $call="cp res_joint_spec.txt res_joint_spec.txt.".$j.".txt";
    system($call);

    #exit;
#    print "$j/$NC $LIST[$j]\n";
}
close(OUT);

$call="ls mod_joint_spec.txt.?.txt > list.now";
system($call);
$call="spectra2img.pl list.now mod_joint_spec.fits";
system($call);
$call="ls model_spec.txt.?.txt > list.now";
system($call);
$call="spectra2img.pl list.now model_spec.fits";
system($call);
$call="ls res_joint_spec.txt.?.txt > list.now";
system($call);
$call="spectra2img.pl list.now res_joint_spec.fits";
system($call);


$call="date > date.end";
system($call);


$out_coeffs="coeffs_".$outfile;
$out_all="all_".$outfile;
$call="joint_loop_table.pl ".$outfile." ".$out_coeffs." ".$out_all;
system($call);



exit;


sub loop {
    my $name=$_[0];
    my $NF=$_[1];
    my $i=0;
    if ($NF<$nf) {
	for ($i=0;$i<$nn[$NF];$i++) {
	    my $name_now=$name." ".$nam_file[$NF][$i];
	    my $val=loop($name_now,$NF+1);
	}
    } else {
	$LIST[$NC]=$name;
	$NC++;
#	print "$name\n";
	return $name;
    }
}

