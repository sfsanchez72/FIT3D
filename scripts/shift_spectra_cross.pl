#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<7) {
    print "USE: shift_spectra_cross.pl INPUT.fits COMPARISON.fits SHIFTED.fits y_index start_index delta_index plot SHIFT_file.txt\n";
    exit;
}

$input_cube=$ARGV[0];
$comparison=$ARGV[1];
$output_cube=$ARGV[2];
$y_index=$ARGV[3];
$w_min=$ARGV[4];
$n_max=$ARGV[5];
$plot=$ARGV[6];
$dist_cor_file=$ARGV[7];


print "Paso $input_cube\n";
$nax=read_naxes($input_cube);   
print "Paso\n";
@naxis=@$nax;
print "$naxis[0] $naxis[1]\n"; 
print "Reading file $input_cube...\n";
#exit;
@in_cube=read_img($input_cube);
print "Readed...\n";

print "Paso $comparison\n";
$nax=read_naxes($comparison);   
print "Paso\n";
@naxis=@$nax;
print "$naxis[0] $naxis[1]\n"; 
print "Reading file $input_cube...\n";
#exit;
@in_comp=read_img($comparison);
print "Readed...\n";


$start_w=0;
$delta_w=1;
$npow=4096;
print "We start to compute the median...\n";
$back_top=0;
$nm=0;
$nm_start=0;

for ($j=0;$j<$naxis[1];$j++) {
    for ($i=0;$i<$naxis[0];$i++) {
	$out_cube[$j][$i]=$in_cube[$j][$i];
    }
}

$min=1e12;
$max=-1e12;


for ($k=0;$k<$naxis[0];$k++) {    
    $w=$start_w+$delta_w*$k;    
    if (($w>$w_min)&&($nm<$n_max)) {
	
	$median[$nm]=$in_comp[$y_index][$k];
	$wave_med[$nm]=$k;
	$wave[$nm]=$nm-$n_max/2;
	if ($nm==0) {
	    $nm_start=$k;
	}
	if ($max<$median[$nm]) {
	    $max=$median[$nm];
	}
	if ($min>$median[$nm]) {
	    $min=$median[$nm];
	}
#	print "$wave[$nm] $median[$nm]\n";
	$nm++;
    }
}
#print "NM=$nm\n";
if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.6);           # Set character height
#    pgenv($w_min,$w_min+$delta_w*$nm,-100,1200,0,0);
#	    pgenv($w_min,$w_min+$delta_w*$nm,$min-1,$max+1,0,0);
#    if ($max<0) {
#	$max=100;
#	$min=-100;
#    }
    pgenv($nm_start+$nm*0.25,$nm_start+0.75*$nm,$min-5,$max+5,0,0);
#	    pgenv($start_w,$start_w+$delta_w*$naxis[2],-100,1200,0,0);
    pgsch(1.6);           # Set character height
    pglabel("X-axis","Counts","Spectral Section");
    pgsch(2.2);           # Set character height
    pgsci(1);
    pgline($nm,\@wave_med,\@median);    
    pgsci(1);
    pgclose;
    pgend;
    print "Press Enter"; <stdin>;
}


print "We start to cross-correlate...\n\n";

#
# We cross-correlate now!
#
$js=0;
$shift_min=1000;
$shift_max=-1000;


$j=$y_index;    
my @spectrum;
$kk=0;
for ($k=$nm_start;$k<($nm_start+$nm);$k++) {    
    $spectrum[$kk]=$in_cube[$j][$k];
    $kk++;
}
#    print "KK=$kk\n";
my $fft1=new Math::FFT(\@spectrum);
if ($js==0) {
#	    if ($reffile eq "") {
#		$fft2=new Math::FFT(\@spectrum);
#	    } else {		
    $kkk=0;
    for ($k=$nm_start;$k<($nm_start+$nm);$k++) {    
	$spectrum_ref[$kkk]=$in_comp[$j][$k];
	$kkk++;
    }
    $fft2=new Math::FFT(\@spectrum_ref);
#	    }
}
#	print "DONE\n";
my $corr= $fft1->correl($fft2);
@a_tmp=@$corr;
#	print "Result\n";
for ($k=0;$k<$nm;$k++) {    
    if ($nm<$n_max/2) {
	$a_corr[$k]=$a_tmp[$k+int($nm/2)];
    } else {
	$a_corr[$k]=$a_tmp[$k-int($nm/2)];
    }
}
$max=-10000;
$min=1e17;
$sum=0;
for ($k=0;$k<$nm;$k++) {    
    if ($a_corr[$k]<$min) {
	$min=$a_corr[$k];
    }
    if ($a_corr[$k]>$max) {
	$max=$a_corr[$k];
	$n_max=$k;
    }
    $sum=$sum+$a_corr[$k];
}
if ($sum>0) {
    $na=$n_max-1;
    $nb=$n_max;
    $nc=$n_max+1;
    if ($na<0) {
	$na=$nm+$na;
    }
    if ($nc>=$nm) {
	$nc=$nc-$nm;
    }
    
    
    $a=$wave[$na];
    $b=$wave[$nb];
    $c=$wave[$nc];
    $fa=-$a_corr[$na];
    $fb=-$a_corr[$nb];
    $fc=-$a_corr[$nc];
#	    $w_max=$b-0.5*((($b-$a)**2)*($fa-$fc)-(($b-$c)**2)*($fb-$fa))/(($b-$a)*($fb-$fc)-($b-$c)*($fb-$fa));
    $den=($fc-2*$fb+$fa);
    if ($den!=0) {
	$w_max=$c-($b-$a)*(($fc-$fb)/$den+0.5);
    } else {
	$w_max=0;
    }
    
    
} else {
    $w_max=0;
}

for ($j=0;$j<$naxis[1];$j++) {    
    if ($plot==1) {
	
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	#pgenv(0,$nm,$min-1,$max+10,0,0);
	pgenv(-$nm/2,$nm/2,$min-0.1,$max+1,0,0);
	pgsch(1.2);           # Set character height
	pglabel("Id.","Counts","[x,y]=($i,$j)");
	pgsch(1.2);           # Set character height
	pgsci(1);
	pgline($nm,\@wave,\@a_corr);
#	    pgsci(3);
#	    pgpoint($nm,\@wave,\@a_corr,2);
	pgsch(2.2);
	pgsci(2);
	pgpoint(3,[$a,$b,$c],[-$fa,-$fb,-$fc],3);
	pgsci(8);
	pgpoint(1,[$w_max],[-$fb],16);
	pgsch(1.2);
	pgsci(1);
	pgclose;
	pgend;
    }
    
    $shift[$js]=$w_max;#+1;
    $index[$js]=$j;
    if ($shift_min>$shift[$js]) {
	$shift_min=$shift[$js];
    }
    
    if ($shift_max<$shift[$js]) {
	$shift_max=$shift[$js];
    }
#	print "$j $js $shift[$js]\n";
    $js++;
    
    
    
}


#
#cleaning
#
my @new_shift=@shift;
my @new_index=@index;
$ks=0;
$nbox=5;
#for ($j=$nbox;$j<$js-$nbox;$j++) {
for ($j=0;$j<$js;$j++) {
    my @tmp;
    $k=0;
    if ($j<$nbox) {
	$j_min=$j;
	$j_max=$j+2*$nbox;
    } else {
	if ($j>=($js-$nbox)) {
	    $j_min=$j-2*$nbox;
	    $j_max=$j;
	} else {
	    $j_min=$j-$nbox;
	    $j_max=$j+$nbox;
	}
    }

#    for ($i=$j-$nbox;$i<$j+$nbox;$i++) {
    for ($i=$j_min;$i<$j_max;$i++) {
	$tmp[$k]=$shift[$i];
	$k++;
    }
    $me=median(@tmp);
    $sig=sigma(@tmp);
#    if (abs($shift[$j]-$me)<$nsigma*$sig) {
	$new_shift[$ks]=$shift[$j];
	$new_index[$ks]=$index[$j];
	$ks++;
#    } else {
#	$shift[$j]=$me;
#    }
}




print "DONE\n";
for ($j=0;$j<$naxis[1];$j++) {
    $index_all[$j]=$j;
    $shift_all[$j]=$shift[$j];
}


#exit;
print "Shifting the new cube...\n";
for ($j=0;$j<$naxis[1];$j++) {
    my @in_spec;
    my @in_wave;
    $min=10000;
    $max=-10000;
#	print "$j $i $shift_all[$j][$i]\n";
    if ($shift_all[$j]!=0) {
	for ($k=0;$k<$naxis[0];$k++) {    
	    $in_spec[$k]=$in_cube[$j][$k];
	    $in_wave[$k]=$k;
#	    $out_wave[$k]=$k-$shift_all[$j];
	    $out_wave[$k]=$k-$shift_all[$j];
	    if ($in_spec[$k]<$min) {
		$min=$in_spec[$k];
	    }
	    if ($in_spec[$k]>$max) {
		$max=$in_spec[$k];
	    }
	    
	}
	my $out_spec_pdl = interpol(pdl(@in_wave), pdl(@out_wave), pdl(@in_spec));
	for ($k=0;$k<$naxis[0];$k++) {
	    $out_spec[$k]=$out_spec_pdl->slice($k)->sclr;		
	    $out_cube[$j][$k]=$out_spec[$k];
	}
	
	
	if ($plot==1) {
	    pgbegin(0,"/xs",1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.6);           # Set character height
#    pgenv($w_min,$w_min+$delta_w*$nm,-100,1200,0,0);
#	    pgenv($w_min,$w_min+$delta_w*$nm,$min-1,$max+1,0,0);
#	    $min=0.5*$max;

	    if ($max<0) {
		$max=100;
		$min=-100;
	    }
	    pgenv($nm_start+$nm*0.25,$nm_start+0.75*$nm,$min-1,$max+500,0,0);
#	    pgenv($start_w,$start_w+$delta_w*$naxis[2],-100,1200,0,0);
	    pgsch(1.6);           # Set character height
	    pglabel("Id.","Counts","($i,$j)");
	    pgsch(2.2);           # Set character height
	    pgsci(1);
	    pgline($nm,\@wave_med,\@median);    
	    pgsci(2);
	    pgline($naxis[0],\@in_wave,\@in_spec);    
	    pgsci(3);
	    pgline($naxis[0],\@out_wave,\@in_spec);    
	    pgsci(5);
	    pgline($naxis[0],\@in_wave,\@out_spec);    
	    pgsci(1);
	    pgclose;
	    pgend;
	    }
    }
    
}

print "DONE\n";
print "Writting the $output_cube\n";
system("rm $output_cube");
system("cp $input_cube $output_cube");
update_fits($output_cube,$nax,2,\@out_cube);
print "DONE\n";
print "Writting the distortion correct at $dist_cor_file";
open(FH,">$dist_cor_file");
for ($j=0;$j<$naxis[1];$j++) {
    print FH "$index_all[$j] $shift_all[$j]\n";
}
close(FH);
print "DONE\n";


exit;





#open(FH,">$output_back");
#$m=1;
#for ($k=0;$k<$npow;$k++) {
#    print FH "$m $w[$k] $back_spec[$k]\n";
#    $m++;
#}
#close(FH);











