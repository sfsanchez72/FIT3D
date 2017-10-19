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


if ($#ARGV<4) {
    print "USE: dist_cor_cross_ascii.pl input_spec.txt ref_spec.txt start_index delta_index plot\n";
    exit;
}

$input_spec=$ARGV[0];
$ref_spec=$ARGV[1];
$w_min=$ARGV[2];
$n_max=$ARGV[3];
$plot=$ARGV[4];
$y_index=0;
$nsigma=2;
if ($#ARGV==6) {
    $y_index=$ARGV[7];
}
if ($#ARGV==7) {
    $y_index=$ARGV[7];
    $nsigma=$ARGV[8];
}
if ($#ARGV==9) {
    $y_index=$ARGV[7];
    $nsigma=$ARGV[8];
}


#print "reading $input_spec\n";
$n1=0;
open(FH,"<$input_spec");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $in_wave[$n1]=$data[0]-1;
    $in_spec[$n1]=$data[2];
    $n1++;
}
#print "Done\n";
close(FH);

#print "reading $ref_spec\n";
$n2=0;
open(FH,"<$ref_spec");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $ref_wave[$n2]=$data[0]-1;
    $ref_spec[$n2]=$data[2];
    $n2++;
}
#print "Done\n";
close(FH);

#if ($n1!=$n2) {
#    print "WArnning N1=$n1!=$n2=N2\n";
#    exit;
#}
$nm=0;
$nm_start=0;
$max=-1e12;
$min=1e12;
for ($k=0;$k<$n1;$k++) {    
    $w=$k;    
    if (($w>$w_min)&&($nm<$n_max)) {
	
	$median[$nm]=$in_spec[$k];
	$median_ref[$nm]=$ref_spec[$k];
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
	if ($max<$median_ref[$nm]) {
	    $max=$median_ref[$nm];
	}
	if ($min>$median_ref[$nm]) {
	    $min=$median_ref[$nm];
	}


#	print "$wave[$nm] $median[$nm]\n";
	$nm++;
    }
}

#$m1=mean(@median);
#$m2=mean(@median_ref);
#for ($k=0;$k<$n1;$k++) {    
#    $ref_spec[$k]=$ref_spec[$k]*($m1/$m2);
#}

#print "NM=$nm\n";
if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.6);           # Set character height
    pgenv($nm_start+$nm*0.25,$nm_start+0.75*$nm,$min-5,$max+5,0,0);
    pgsch(1.6);           # Set character height
    pglabel("X-axis","Counts","Spectral Section");
    pgsch(2.2);           # Set character height
    pgsci(1);
    pgline($nm,\@wave_med,\@median);    
    pgsls(2);
    pgsci(2);
    pgline($n1,\@in_wave,\@in_spec);    
    pgsls(1);
    pgsci(3);
    pgline($n2,\@ref_wave,\@ref_spec);    
    pgsci(1);
    pgclose;
    pgend;
    print "Press Enter"; <stdin>;
}
#exit;

#print "We start to cross-correlate...\n\n";

#print "NSTART=$nm_start\n";
#
# We cross-correlate now!
#
$js=0;
$shift_min=1000;
$shift_max=-1000;
my @spectrum;
my @spectrum_ref;
$kk=0;
for ($k=$nm_start;$k<($nm_start+$nm);$k++) {    
    $spectrum[$kk]=$in_spec[$k];
    $spectrum_ref[$kk]=$ref_spec[$k];
    $kk++;
}
my $fft1=new Math::FFT(\@spectrum);
$fft2=new Math::FFT(\@spectrum_ref);
my $corr= $fft1->correl($fft2);
@a_tmp=@$corr;
#print "N=$nm, $n_max\n";
for ($k=0;$k<$nm;$k++) {    
    if ($nm<($n_max/2)) {
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
    $den=($fc-2*$fb+$fa);
    if ($den!=0) {
	$w_max=$c-($b-$a)*(($fc-$fb)/$den+0.5);
    } else {
	$w_max=0;
    }
    
	
} else {
    $w_max=0;
}
    
#$w_max=$w_max-$n_max/2;
    
if ($plot==1) {
    
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    #pgenv(0,$nm,$min-1,$max+10,0,0);
    pgenv(-$nm/2,$nm/2,$min-0.1,$max+1,0,0);
    pgsch(1.2);           # Set character height
    pglabel("Id.","Counts","");
    pgsch(1.2);           # Set character height
    pgsci(1);
    pgline($nm,\@wave,\@a_corr);
#    pgsci(2);
#    pgline($nm,\@in_wave,\@a_tmp);
#	    pgsci(3);
#	    pgpoint($nm,\@in_wave,\@a_corr,2);
    pgsch(2.2);
    pgsci(2);
    pgpoint(3,[$a,$b,$c],[-$fa,-$fb,-$fc],3);
    pgsci(8);
    pgpoint(1,[$w_max],[-$fb],16);
    pgsch(1.2);
    pgsci(1);
	pgclose;
    pgend;
    print "Press Enter"; <stdin>;
}

$max=-1e12;
$min=1e12;

for ($k=0;$k<$n1;$k++) {    
    $out_wave[$k]=$k-$w_max;
    if ($in_spec[$k]<$min) {
	$min=$in_spec[$k];
    }
    if ($in_spec[$k]>$max) {
	$max=$in_spec[$k];
    }
}
my $out_spec_pdl = interpol(pdl(@in_wave), pdl(@out_wave), pdl(@in_spec));
for ($k=0;$k<$n1;$k++) {
    $out_spec[$k]=$out_spec_pdl->slice($k)->sclr;		
    if ($max<$out_spec[$k]) {
	$max=$out_spec[$k];
    }
    if ($min>$out_spec[$k]) {
	$min=$out_spec[$k];
    }
#    print "$k $in_wave[$k] $in_spec[$k] $out_spec[$k]\n";
}
	

if ($max<0) {
    $max=100;
    $min=-100;
}

if ($plot==1) {    
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv($nm_start+$nm*0.25,$nm_start+0.75*$nm,$min-1,$max+10,0,0);
    pgsch(1.6);           # Set character height
    pglabel("Id.","Counts","");
    pgsch(2.2);           # Set character height
    pgsci(1);
    pgline($nm,\@wave_med,\@median);    
    pgsci(2);
    pgline($n1,\@in_wave,\@in_spec);    
    pgsci(3);
    pgline($n1,\@out_wave,\@in_spec);    
    pgsci(5);
#    pgsci(3);
    pgline($n2,\@ref_wave,\@ref_spec);    
#    pgline($n1,\@in_wave,\@out_spec);    
    pgsci(1);
    pgclose;
    pgend;
    print "Press Enter"; <stdin>;
}

print "$input_spec $w_max\n";

exit;





