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
    print "USE: auto_ssp_elines_several_Av_log_cube.pl SPEC1.RSS.fits BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] \n";
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

$outfile=$ARGV[2];
$h=$pdl->gethdr;
$mod_cube_file="mod_".$outfile.".fits";
$res_cube_file="res_".$outfile.".fits";
$gas_cube_file="gas_".$outfile.".fits";
$no_gas_cube_file="no_gas_".$outfile.".fits";

$mod_cube=zeroes($nx,$ny,$nz);
$res_cube=zeroes($nx,$ny,$nz);
$gas_cube=zeroes($nx,$ny,$nz);
$no_gas_cube=zeroes($nx,$ny,$nz);

$mod_cube->sethdr( $h );
$res_cube->sethdr( $h );
$gas_cube->sethdr( $h );
$no_gas_cube->sethdr( $h );






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
open(FL,"<F.limit");
$FL=<FL>;
chop($FL);
close(FL);

for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$nz_med=int($nz/2);
	$F=$pdl->at($i,$j,$nz_med);
	if (($F>$FL)&&($F<1e30)) {
	    $call="rm -f spec.auto.txt";
	    system($call);
	    $call="rm -f tmp.out";
	    system($call);
	    $call="get_spec_cube.pl ".$infile." ".$i." ".$j." spec.auto.txt";
	    system($call);
#	print "$call\n Press Enter\n"; <stdin>;
	    $call="auto_ssp_elines_several_Av_log.pl spec.auto.txt ";
	    for ($jj=1;$jj<$NP;$jj++) {
		$call=$call." ".$PARAM[$jj];
	    }
	    $call=$call." >> auto.log";
	    system($call);
#	    print "$call\n Press Enter\n"; <stdin>;
	    open(FH,"<tmp.out");
	    $line=<FH>;
	    if ($i==0) {
		print OUT "$line";
	    }
	    $line=<FH>;
	    chop($line);
	    print OUT "$nx $ny $i $j $line\n";
	    close(FH);
	    open(ELINES,">>$elines");
	    $ne=0;
	    open(IN,"<elines_tmp.out");
	    while($in=<IN>) {
		if ($in =~ "eline") {
		    $line_in[$ne]=$in;
		    $ne++;
		}
	    }
	    close(IN);
	    print ELINES "$nx $ny $i $j\n";
	    print ELINES "$ne\n";
	    for ($ii=0;$ii<$ne;$ii++) {
		print ELINES "$line_in[$ii]";
	    }
	    close(ELINES);
#	$call="cat elines_tmp.out >> ".$elines;
#	system($call);
	    $call="cat coeffs.out >> ".$coeffs;
	    system($call);



	    add_file_mod($i,$j,"org_spec.txt");
	    add_file_res($i,$j,"res_joint_spec.txt");
	    add_file_gas($i,$j,"gas_spec.txt");
	    add_file_no_gas($i,$j,"no_gas_spec.txt");
	    





	print "$i/$nx $j/$ny DONE\n";
	} else {
	print "$i/$nx $j/$ny PASSED BY\n";
	}
    
    }
}
close(OUT);

$mod_cube->wfits($mod_cube_file);
$res_cube->wfits($res_cube_file);
$gas_cube->wfits($gas_cube_file);
$no_gas_cube->wfits($no_gas_cube_file);



exit;



sub add_file_mod {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($mod_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}




sub add_file_res {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($res_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}

sub add_file_gas {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($gas_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}

sub add_file_no_gas {
    my $I=$_[0];
    my $J=$_[1];
    my $file=$_[2];
    my $i,$line;
    $i=0;
    open(FH,"<$file");
    while($line=<FH>) {
        my @data=split(" ",$line);
        my $val=$data[2];
        set($no_gas_cube,$I,$J,$i,$val);
        $i++;
    }
    close(FH);
    return;
}






