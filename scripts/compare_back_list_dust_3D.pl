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



$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<12) {
    print "USE: compare_back_list_dust_3D.pl RSS.fits BACK_LIST OUTFILE_FULL OUTFILE_AGEMET MASK_LIST REDSHIFT SIGMA_DISP WAVE_SCALE FLUX_THRESHOLD PLOT Av_min Av_max delta_Av  [min max] [wmin wmax]\n";
    exit;
}

$unc_file=$ARGV[0];
$back_list=$ARGV[1];
$outfile=$ARGV[2];
$out_AM=$ARGV[3];
$mask_list=$ARGV[4];
$redshift=$ARGV[5];
$sigma=$ARGV[6];
$wave_scale=$ARGV[7];
$out_file="junk.junk";
$factor=1;
$box=1;
$flux_limit=$ARGV[8];
$plot=$ARGV[9];
$Av_min=$ARGV[10];
$Av_max=$ARGV[11];
$delta_Av=$ARGV[12];



$smooth=1;

$def=0;
if ($#ARGV==14) {
    $min=$ARGV[13];
    $max=$ARGV[14];
    $def=1;
}

if ($#ARGV==16) {
    $min=$ARGV[13];
    $max=$ARGV[14];
    $min_wave=$ARGV[15];
    $max_wave=$ARGV[16];
    $def=2;
}

$pdl=rfits($unc_file);
($nx,$ny)=$pdl->dims;
$res_pdl=$pdl;
$mod_pdl=zeroes($nx,$ny);


open (OUT,">$outfile"); 
open (AM,">$out_AM"); 
for ($i=0;$i<$ny;$i++) {
    print "Processing $i/$ny...";
    $call="img2spec.pl ".$unc_file." ".$i." spec_junk.txt";
    system($call);
    print OUT "$i ";
    $call="elines_find.pl spec_junk.txt 4 0.05 elines_find.junk";
    system($call);
    $zmin=0;
    $zmax=$redshift+0.02;
    $call="redshift_find.pl elines_find.junk emission_lines.detected 3 ".$zmin." ".$zmax;
    system($call);
    open(RED,"<redshift_find.out");
    $RED=<FH>;
    chop($RED);
    ($red_now,$ered_now,$match)=split(" ",$RED);
    close(RED);
    if ($red_now==0) {
	$red_now=$redshift;
	$ered_now=0;
    }

#
# We check the intensity level
#
    $call="spec_plot.pl spec_junk.txt /xs > spec_plot.out";
    system($call);
    open(FLUX,"<spec_plot.out");
    $line=<FLUX>;
    chop($line);
    @data=split(" ",$line);
    $median_flux=$data[0];
    close(FLUX);

#    print "PASO\n";
    print OUT "$median_flux $flux_limit\n";
#    print "PASO $median_flux $flux_limit\n";
    if ($median_flux>$flux_limit) {
	
	$call="compare_back_list_dust.pl spec_junk.txt ".$back_list." comp_back_junk.out ".$mask_list." ".$red_now." ".$sigma." ".$wave_scale." ".$plot." ".$Av_min." ".$Av_max." ".$delta_Av;
	#print "$call\n";
	if ($def==1) {
	    $call=$call." ".$min." ".$max;
	}
	if ($def==2) {
	    $call=$call." ".$min." ".$max." ".$min_wave." ".$max_wave;
	}
	system($call);
	
	open(FH,"<median.flux");
	$median_flux=<FH>;
	chop($median_flux);
	#print "MEDIAN_FLUX=$median_flux ";
	close(FH);


#    system($call);
#    print "\n $call\n";
 #
# We copy the result
#
	open(FH,"<comp_back_junk.out");
	while($line=<FH>) {
	    chop($line);
	    print OUT "$line\n";
	}
	close(FH);
	
#	$call="plot_Age_Met.pl comp_back_junk.out ";
#	if ($plot==1) {
#	    $call=$call." /xs";
#	} else {
#	    $call=$call." /null";
#	}

#
# We copy the residual and model spectra
#
	$n=0;
	open(FH,"<res_spec.txt");
	while ($line=<FH>) {
	    chop($line);
	    @data=split(" ",$line);	    
	    set($res_pdl,($n,$i),$data[2]);
	    $n++;
	}

	$n=0;
	open(FH,"<model_spec.txt");
	while ($line=<FH>) {
	    chop($line);
	    @data=split(" ",$line);
	    set($mod_pdl,($n,$i),$data[2]);
	    $n++;
	}

    }



    open(FH,"<compare_back_list_dust.out");
    $line=<FH>;
    chop($line);
    if ($median_flux>$flux_limit) {
	@data=split(" ",$line);
	print AM "$median_flux AGEMET $data[1] 0.0 $data[2] 0.0 $data[3] 0.0 $red_now $ered_now\n";
    } else {
	print AM "$median_flux AGEMET 0.0 0.0 0.0 0.0 $red_now $ered_now\n";
    }
    close(FH);
    print "DONE\n";
}
close(AM);
close(OUT);

$h=$pdl->gethdr;
$res_pdl->wfits("res_spec.fits");
$mod_pdl->sethdr($h);
$mod_pdl->wfits("model_spec.fits");




exit;
