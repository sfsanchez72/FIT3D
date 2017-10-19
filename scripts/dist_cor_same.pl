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


if ($#ARGV<2) {
    print "USE: dist_cor_same.pl EXTRACTED.fits CORRECTED.fits delta_index\n";
    exit;
}

$input_cube=$ARGV[0];
$output_cube=$ARGV[1];
$output_ps=$ARGV[1];
#$dist_cor_file=$ARGV[2];
$output_ps =~ s/.fits/.ps/g;
#$smooth=$ARGV[3];
#$w_min=$ARGV[3];
$delta_n=$ARGV[2];
$plot=0;


print "Reading $input_cube\n";
$nax=read_naxes($input_cube);   
#print "DINE\n";
@naxis=@$nax;
#print "$naxis[0] $naxis[1]\n"; 
print "Reading file $input_cube...\n";
#exit;
@in_cube=read_img($input_cube);
print "DONE...\n";

print "Reading the distortion correction\n";
$n=0;
$shift_min=1e12;
$shift_max=-1e12;
for ($j=0;$j<$naxis[1];$j++) {
    $index_all[$j]=$j;
    $shift_all[$j]=$delta_n;
    if ($shift_min>$shift_all[$j]) {
	$shift_min=$shift_all[$j];
    }
    if ($shift_max<$shift_all[$j]) {
	$shift_max=$shift_all[$j];
    }
}
$n=$naxis[1];
print "N=$n\n";

if ($n!=$naxis[1]) {
    print "The dimension of the distortion correction file, $dist_cor_file,\n";
    print "and the extracted image, $input_cube, do not match\n";
    exit;
}

#($start_w,$delta_w)=read_img_headers($input_cube,["CRVAL3","CDELT3"]);
$start_w=0;
$delta_w=1;
$npow=4096;
#print "$start_w,$delta_w\n";
#exit;
print "We start to compute the median...\n";
$back_top=0;
$nm=0;
$nm_start=0;

for ($j=0;$j<$naxis[1];$j++) {
    for ($i=0;$i<$naxis[0];$i++) {
	$out_cube[$j][$i]=$in_cube[$j][$i];
    }
}

#$w=$start_w+$delta_w*$k;
$min=1e12;
$max=-1e12;

for ($k=0;$k<$naxis[0];$k++) {    
    $w=$start_w+$delta_w*$k;    
    if (($w>$w_min)&&($nm<$n_max)) {
	$median[$nm]=$in_cube[2][$k];
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
print "NM=$nm\n";
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


if ($plot==1) {
    
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.6);           # Set character height
    pgenv(0,$naxis[1],$shift_min-1,$shift_max+1,0,0);
    pgsch(1.6);           # Set character height
    pglabel("Y index","Delta","Shifts");
    pgsch(2.2);           # Set character height
    pgsci(1);
#    pgpoint($js,\@index,\@shift,1);
#     pgsci(3);
#    pgpoint($ks,\@new_index,\@new_shift,22);
    pgsci(2);
    pgline($naxis[1],\@index_all,\@shift_all);
    pgsci(1);
    pgclose;
    pgend;
    print "Press Enter"; <stdin>;    
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
    print "J=$j/$n\n";
}

print "DONE\n";
print "Writting the $output_cube\n";
system("rm $output_cube");
system("cp $input_cube $output_cube");
update_fits($output_cube,$nax,2,\@out_cube);
print "DONE\n";

exit;





#open(FH,">$output_back");
#$m=1;
#for ($k=0;$k<$npow;$k++) {
#    print FH "$m $w[$k] $back_spec[$k]\n";
#    $m++;
#}
#close(FH);











