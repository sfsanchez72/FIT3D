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
    print "USE: multi_SSP.pl SPEC1.txt MULTI_SSP.list OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
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
for ($j=0;$j<$#ARGV+1;$j++) {
    $PARAM[$j]=$ARGV[$j];
}
$NP=$#ARGV+1;
$PARAM[1]="tmp.ssp.fits";
$PARAM[2]="tmp.out";
$multi_ssp=$ARGV[1];
#
# We establish the number of templates
#

$nt=0;
$NT=1;
open(FH,"<$multi_ssp");
while($line=<FH>) {
    chop($line);
    $list[$nt]=$line;
    $nwc=0;
    open(TP,"<$line");
    while($tp=<TP>) {
	chop($tp);
	$SSP[$nt][$nwc]=$tp;
	#yprint "$nt $nwc $SSP[$nt][$nwc]\n";
	for ($i=0;$i<$nt;$i++) {
	    $M[$i][$nwc]=$M[$i][$nwc]."$nwc";
	}
	$nwc++;
    }
    close(TP);
    $ns[$nt]=$nwc;   
    $NT=$NT*$nwc;
    $nt++;
}
close(FH);
($ns_min,$ns_max)=minmax(@ns);

$add="";
$add=COMBI(0,$add);
print "$add\n";





#
# COMBINATORIA....
#

exit;


sub COMBI {
    my $it=$_[0];
    my $add=$_[2];
    my $i,$j;
    if ($it<$nt) {
	for ($j=0;$j<$ns[$it];$j++) {    
	    my $oadd=$add." ".$SSP[$it][$j];
	    $add=COMBI($it+1,$oadd);
	    print "$it/$nt $j $add\n";

	}
    }
    return $add;	    
}


$elines="elines_".$ARGV[2];
$call="rm -f ".$elines;
$coeffs="coeffs_".$ARGV[2];
$call="rm -f ".$elines;
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
    open(FH,"<spec.auto.txt");
    open(TMP,">spec.auto.tmp");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$ef=0.1*$data[2];
	print TMP "$line $ef 1\n";
    }
    if ($no_phot==1) {
	for ($j=0;$j<$np;$j++) {
	    $wp=$pdl_phot->at($j,$nyp);
	    $fp=$pdl_phot->at($j,$i);
	    $efp=0.05*$fp;
	    print TMP "$j $wp $fp $efp 2\n";
	}
    }


    close(TMP);
    close(FH);

#    $call="auto_ssp_elines_several_Av_log.pl spec.auto.txt ";
    $call="auto_ssp_elines_several_Av_log.pl spec.auto.tmp ";
    for ($j=1;$j<$NP;$j++) {
	$call=$call." ".$PARAM[$j];
    }
#    print "$call\n"; exit;
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


    add_file($i,0,"org_spec.txt");
    add_file($i,1,"model_spec.txt");
    add_file($i,2,"mod_joint_spec.txt");
    add_file($i,3,"res_spec.txt");
    add_file($i,4,"res_joint_spec.txt");
    add_file($i,5,"no_gas_spec.txt");
}
close(OUT);

$pdl_output->wfits($output_file);


exit;


sub add_file {
    my $J=$_[0];
    my $K=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
	my @data=split(" ",$line);
	my $val=$data[2];
	set($pdl_output,$i,$J,$K,$val);
	$i++;
    }
    close(FH);
    return;
}
