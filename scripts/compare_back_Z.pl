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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<9) {
    print "USE: compare_back_Z.pl SPEC1.txt BACK_SPEC OUTFILE MASK_LIST VEL_MIN VEL_MAX DELTA_VEL SIGMA_DISP WAVE_SCALE PLOT [min max] [wmin wmax]\n";
    exit;
}

$unc_file=$ARGV[0];
$c_file=$ARGV[1];
$outfile=$ARGV[2];
$mask_list=$ARGV[3];
$vel_min=$ARGV[4];
$vel_max=$ARGV[5];
$delta_vel=$ARGV[6];
$sigma=$ARGV[7];
$wave_scale=$ARGV[8];
$out_file="junk.junk";
$factor=1;
$box=1;
$plot=$ARGV[9];
$smooth=1;
$MIN_CHISQ=1e12;

$def=0;
if ($#ARGV==11) {
    $min=$ARGV[10];
    $max=$ARGV[11];
    $def=1;
}

if ($#ARGV==13) {
    $min=$ARGV[10];
    $max=$ARGV[11];
    $min_wave=$ARGV[12];
    $max_wave=$ARGV[13];
    $def=2;
}




if ($mask_list eq "none") {
    $nmask=0;
} else {
    open(FH,"<$mask_list");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$start_mask[$nmask]=$data[0];
	$end_mask[$nmask]=$data[1];
	$nmask++;
    }
    close(FH);
}




$n_unc=0;
$y_min=1e12;
$y_max=-1e12;
open(FH,"<$unc_file");
$i_scale=0;
while($line=<FH>) {
    @data=split(" ",$line);
    $index_unc[$n_unc]=$data[0];
    $wave_unc[$n_unc]=$data[1];
    $flux_unc[$n_unc]=$data[2];    
    if ($flux_unc[$n_unc]<$y_min) {
	$y_min=$flux_unc[$n_unc];
    }
    if ($flux_unc[$n_unc]>$y_max) {
	$y_max=$flux_unc[$n_unc];
    }
    if ($n_unc>0) {
	if (($wave_unc[$n_unc-1]<$wave_scale)&&($wave_unc[$n_unc]>$wave_scale)) {
	    $i_scale=$n_unc;	    
	}
    }
    $masked[$n_unc]=1;
    $w_test=$wave_unc[$n_unc-1];
    for ($j=0;$j<$nmask;$j++) {
	if (($w_test>$start_mask[$j])&&($w_test<$end_mask[$j])) {
	    $masked[$n_unc]=0;	    
	}
    }
    
    if ($def==2) {
	if ($w_test<$min_wave) {
	    $masked[$n_unc]=0;	    
	}
	if ($w_test>$max_wave) {
	    $masked[$n_unc]=0;	    
	}
    }

    $flux_masked[$n_unc]=$flux_unc[$n_unc]*$masked[$n_unc];

#    if ($def==2) {
#	if (($data[1]>$min_wave)&&($data[1]<$max_wave)) {
#	    $n_unc++;
#	}
#    } else {
	$n_unc++;
#    }
}
close(FH);
if ($def==2) {
    $y_min=$min;
    $y_max=$max;
} else {
    $min_wave=$wave_unc[0];
    $max_wave=$wave_unc[$n_unc-1];
}

if ($def==1) {
    $y_min=$min;
    $y_max=$max;
}
#
# We check the median value
#
$median=median(@flux_masked);
open(MOUT,">median.flux");
print MOUT "$median\n";
close(MOUT);


$dpix_unc=$wave_unc[1]-$wave_unc[0];

open(OUT,">$outfile");

$vel_now=$vel_min;
$chi_sq_min=1e12;
while ($vel_now<$vel_max) {
#for ($iii=0;$iii<$nf;$iii++) {
#    $c_file=$file[$iii];
#    print "$c_file\n";

    $n_c=0;
    open(FH,"<$c_file");
    while($line=<FH>) {
	@data=split(" ",$line);
	$redshift=$vel_now/300000;
	#if (($data[1]>$min_wave)&&($data[1]<$max_wave)) {
	    $wave_c[$n_c]=$data[1]*(1+$redshift);
	    $flux_c_ini[$n_c]=$data[2];
	    $n_c++;
	#}
    }
    $dpix_c=$wave_c[1]-$wave_c[0];
    close(FH);

    
#
# We convolve with a gaussian
# 
#$sigma=10;
    $rsigma=$sigma/$dpix_c;
    @flux_c=@flux_c_ini;
    $box=int(3*$rsigma);
    if ($box>3) {
	for ($i=$box;$i<($n_c-$box);$i++) {
	    $norm=0;
	    $flux_c[$i]=0;
	    for ($j=$i-$box;$j<($i+$box);$j++) {
		$gaus=exp(-0.5*((($i-$j)/$rsigma)**2));
	    $flux_c[$i]=$flux_c[$i]+$flux_c_ini[$j]*$gaus;
		$norm=$norm+$gaus;
	    }
	    $flux_c[$i]=$flux_c[$i]/$norm;
	}
    }
    
    
#print "$n_unc $n_c\n";
#
# We interpolate
#
    my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), pdl(@flux_c));

    $scale=$flux_unc[$i_scale]/($out_spec_pdl->at($i_scale));
    $out_spec_pdl=$out_spec_pdl*$scale;
    
    $min=1e12;
    $max=-1e12;
    $min_rat=1e30;
    $max_rat=-1e30;
    $chi=0;
    for ($j=0;$j<$n_unc;$j++) {
	$out_spec[$j]=$out_spec_pdl->at($j);		
#	print "$j $flux_unc[$j] $chi\n";
	if ($flux_unc[$j]!=0) {
	    if (($wave_unc[$j]>$min_wave)&&($wave_unc[$j]<$max_wave)) {
		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($flux_unc[$j]);
	    }
	}
    }
    $chi_sq=(($chi/$n_unc)**0.5);
    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
	for ($j=0;$j<$n_unc;$j++) {
	    $out_spec[$j]=$out_spec_pdl->at($j);		
	    $res_spec[$j]=$flux_unc[$j]-$out_spec[$j];
	    $model_spec[$j]=$out_spec[$j];
	}
	
    }

#    $chi_sq=($chi)**0.5;
    @data=split(/\//,$c_file);
    $name=$data[$#data];
    
#    print "$name $chi_sq $scale\n";
    print OUT "$name $chi_sq $vel_now\n";

    if ($chi_sq_min>$chi_sq) {
	$velocity=$vel_now;
	$chi_sq_min=$chi_sq;

    }
    
    if ($plot==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	#pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	pgenv($min_wave,$max_wave,$y_min,$y_max,0,0);
	pgsch(1.2);           # Set character height
	pglabel("Wavelength","Counts","");
	pgsci(1);
	pgline($n_unc,\@wave_unc,\@flux_unc);    
	pgsci(2);
	pgline($n_unc,\@wave_unc,\@out_spec);    
	pgsci(3);
	pgline($n_unc,\@wave_unc,\@flux_masked);    
	pgsci(1);
	pgclose;
	pgend;
#    print "Press Enter"; <stdin>;
    }
 $vel_now=$vel_now+$delta_vel;   
}
close(OUT);

open(OUT,">COMPARE_BACK_Z.out");
print OUT "$velocity $chi_sq_min\n";
close(OUT);
print "$velocity $chi_sq_min\n";



exit;

