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
    print "USE: auto_ssp_elines_several_Av_log_rss_new.pl SPEC1.RSS.fits BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
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
$error_infile="e_".$ARGV[0];
$pdl=rfits($infile);
$pdl_error=rfits($infile);
($nx,$ny)=$pdl->dims();
$crval=$pdl->hdr->{CRVAL1};
$crpix=$pdl->hdr->{CRPIX1};
$cdelt=$pdl->hdr->{CDELT1};
for ($j=0;$j<$#ARGV+1;$j++) {
    $PARAM[$j]=$ARGV[$j];
}
$NP=$#ARGV+1;
$PARAM[2]=" tmp.out";
$elines="elines_".$ARGV[2];
$call="rm -f ".$elines;
system($call);
$coeffs="coeffs_".$ARGV[2];
$all="all_".$ARGV[2];
$call="rm -f ".$coeffs;
system($call);
$output_file="output.".$ARGV[2].".rss.fits";
$call="rm -f ".$output_file;
system($call);
$call="rm -f auto.log";
system($call);
open(OUT,">$ARGV[2]");
open(OUT2,">$all");
my $pdl_output=zeroes($nx,$ny,6);
$h = {NAXIS=>3, NAXIS1=>$nx, NAXIS2=>$ny, NAXIS3=>6 , COMMENT=>"output FITS file",
      NAME0 => "org_spec",
      NAME1 => "model_spec",
      NAME2 => "mod_joint_spec",
      NAME3 => "gas_spec",
      NAME4 => "res_joint_spec",
      NAME5 => "no_gas_spec",
      CRVAL3 => $crval,
      CRPIX3 => $crpix,
      CDELT3 => $cdelt
      }; 
$$h{FILENAME} =$output_file;
$pdl_output->sethdr( $h );

#
# Photometric points to add?
#
$phot_file=$pdl->hdr->{PHOTFILE};
$no_phot=0;
if ($phot_file ne "") {
    $no_phot=1;
    $pdl_phot=rfits($phot_file);
    ($np,$ndim)=$pdl_phot->dims;
    $nyp=$ndim-1;
    if ($nyp!=$ny) {
	print "PHOT and SPEC dimensions do not match ($nyp!=$ny)\n"; exit;
    }
}


for ($i=0;$i<$ny;$i++) {
    $call="rm -f spec.auto.txt";
    system($call);
    $call="rm -f tmp.out";
    system($call);
    $call="img2spec_e.pl ".$infile." ".$i." spec.auto.txt";
    system($call);
    open(FH,"<spec.auto.txt");
    open(TMP,">spec.auto.tmp");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$ef=0.1*abs($data[2]);
	if ($ef==0) {
	    $ef=1;
	}
#	print TMP "$line $ef 1\n";
	print TMP "$line 1\n";
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

# NO PHOT
#    $call="auto_ssp_elines_several_Av_log.pl spec.auto.txt ";
    $call="auto_ssp_elines_several_Av_log_new.pl spec.auto.tmp ";
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

    open(EE,">>$elines");
    print EE "#ID $i\n";
    close(EE);
    $call="cat elines_tmp.out >> ".$elines;
    system($call);
    $call="cat coeffs.out >> ".$coeffs;
    system($call);

    open(FH,"<all_tmp.out");
    while ($line=<FH>) {
	chop($line);
	if ($line !~ "#") {
	    print OUT2 "$line";
	} else {
	    if ($i==0) {
		print OUT2 "$line";
	    }
	}
    }


    close(FH);

    print "$i/$ny DONE\n";


    add_file($i,0,"org_spec.txt");
    add_file($i,1,"model_spec.txt");
    add_file($i,2,"mod_joint_spec.txt");
    add_file($i,3,"gas_spec.txt");
    add_file($i,4,"res_joint_spec.txt");
    add_file($i,5,"no_gas_spec.txt");
}
close(OUT2);
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
