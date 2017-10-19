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


#use POSIX;

$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");



if ($#ARGV<4) {
    print "USE: check_eline.pl EXTRACTED.fits start_index delta_index nsearch plot \n";
    exit;
}

$input_cube=$ARGV[0];
$w_guess=$ARGV[1];
$delta_n=$ARGV[2];
$nsearch=$ARGV[3];
$plot=$ARGV[4];

#print "Paso $input_cube\n";
$nax=read_naxes($input_cube);   
#print "Paso\n";
@naxis=@$nax;
#print "$naxis[0] $naxis[1]\n"; 
#print "Reading file $input_cube...\n";
#exit;
@in_cube=read_img($input_cube);

#print "Readed...\n";
$start_w=0;
$delta_w=1;
$npow=4096;
#print "We start to compute the median...\n";
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
$jmin=$w_guess-$delta_n;
$jmax=$w_guess+$delta_n;
for ($nm=0;$nm<$naxis[0];$nm++) {    
    $spec[$nm]=$in_cube[$n_fib][$nm];
    $wave[$nm]=$nm;
    if (($nm>=$jmin)&&($nm<=$jmax)) {
	if ($max<$spec[$nm]) {
	    $max=$spec[$nm];
	}
	if ($min>$spec[$nm]) {
	    $min=$spec[$nm];
	}
    }
#	print "$wave[$nm] $median[$nm]\n";
}
#print "We start look for the 1st peak ...\n\n";


$j=int($jmin);
$peak=0;
$max_peak=-1e12;
$i_peak[$n_fib]=$w_guess;
$w_peak_sum=0;
$spec_sum=0;
while ($j<$jmax) {
    if ($nsearch!=0) {
	$peak=1;
	for ($i=0;$i<$nsearch;$i++) {
	    if ($spec[$j-$i]<$spec[$j-$i-1]) {
		$peak=0;
	    }
	    if ($spec[$j+$i]<$spec[$j+$i+1]) {
	    $peak=0;
	    }
	    
	} 
	if (($peak==1)&&($max_peak<$spec[$j])) {
	    $i_peak[$n_fib]=$j;    
	    $max_peak=$spec[$j];
	}		



	#print "'$j' $w_peak_sum $spec_sum $spec[$j]\n";


    } else {
	if ($max_peak<$spec[$j]) {
	    $i_peak[$n_fib]=$j;    
	    $max_peak=$spec[$j];
	}
    }
    $j++;
}

$na=$i_peak[$n_fib]-1;
$nb=$i_peak[$n_fib];
$nc=$i_peak[$n_fib]+1;
$a=$wave[$na];
$b=$wave[$nb];
$c=$wave[$nc];
$fa=-$spec[$na];
$fb=-$spec[$nb];
$fc=-$spec[$nc];
$den=($fc-2*$fb+$fa);
if ($spec_sum>0) {
    $w_peak_mean[$n_fib]=$w_peak_sum/$spec_sum;
} else {
    $w_peak_mean[$n_fib]=$w_peak_sum;
}
if ($den!=0) {
    $w_peak[$n_fib]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
} else {
    $w_peak[$n_fib]=$w_guess;
}
$flux_peak[$n_fib]=$spec[$nb];

if ($center_index==0) {
    $center_index=$w_peak[$n_fib];
}

#print "AAA $i_peak[$n_fib] '$w_peak[$n_fib]' $w_peak_mean[$n_fib]\n";

#print "NM=$nm\n";
if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.6);           # Set character height
    pgenv($jmin,$jmax,$min-5,$max+5,0,0);
    pgsch(1.6);           # Set character height
    pglabel("X-axis","Counts","Spectral Section");
    pgsch(2.2);           # Set character height
    pgsci(1);
    pgline($naxis[0],\@wave,\@spec);    
    pgsci(7);
#    print "$i_peak[$n_fib]\n";
    $x=$i_peak[$n_fib];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],1);
    pgsci(4);
    $x=$w_peak[$n_fib];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],2);
    pgsci(8);
    $x=$w_peak_mean[$n_fib];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],2);
#    pgptxt();
    pgsci(1);

    pgclose;
    pgend;
#    print "Press Enter"; <stdin>;
	if (($command ne "A")&&($command ne "0")) {
	    print "Press Enter (A=Automatic, 0=No plot)";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }
	    if ($command eq "N") {
		$plot=2;
	    }
	}

}

@Spec=@spec;
$W_guess=$w_guess;

#print "We look for the peaks in all the spectra\n";
$Jmin=$jmin;
$Jmax=$jmax;
$min=1e12;
$max=-1e12;
for ($k=($n_fib+1);$k<$naxis[1];$k++) {
    $w_guess=$w_peak[$k-1];
    $jmin=$W_guess-$delta_n;
    $jmax=$W_guess+$delta_n;

    for ($nm=0;$nm<$naxis[0];$nm++) {    
	$spec[$nm]=$in_cube[$k][$nm];
	$wave[$nm]=$nm;
	if (($nm>=$jmin)&&($nm<=$jmax)) {
	    if ($max<$spec[$nm]) {
		$max=$spec[$nm];
	    }
	    if ($min>$spec[$nm]) {
		$min=$spec[$nm];
	    }
	}
    }


    $peak=0;
    $max_peak=-1e12;
    $j=$jmin;
    $i_peak[$k]=$w_guess;
    $w_peak_sum=0;
    $spec_sum=0;
    while ($j<$jmax) {
	if ($nsearch!=0) {
	    $peak=1;
	    for ($i=0;$i<$nsearch;$i++) {
		if ($spec[$j-$i]<$spec[$j-$i-1]) {
		    $peak=0;
		}
		if ($spec[$j+$i]<$spec[$j+$i+1]) {
		    $peak=0;
		}
		$w_peak_sum=$w_peak_sum+$spec[$j-$i]*($j-$i);
		$spec_sum=$spec_sum+$spec[$j-$i];
#		$SS=$spec[$j-$i];
#		print "$jmin $j $i $SS $w_peak_sum=$w_peak_sum+$spec[$j-$i]*($j-$i) $spec_sum=$spec_sum+$spec[$j-$i]\n";


	    } 
	    if (($peak==1)&&($max_peak<$spec[$j])) {
		$i_peak[$k]=$j;    
		$max_peak=$spec[$j];
	    }
	} else {
	    if ($max_peak<$spec[$j]) {
		$i_peak[$k]=$j;    
		$max_peak=$spec[$j];
	    }
	}
	$j++;	
    }

$na=$i_peak[$k]-1;
$nb=$i_peak[$k];
$nc=$i_peak[$k]+1;
$a=$wave[$na];
$b=$wave[$nb];
$c=$wave[$nc];
$fa=-$spec[$na];
$fb=-$spec[$nb];
$fc=-$spec[$nc];
$den=($fc-2*$fb+$fa);

    if ($spec_sum>0) {
	$w_peak_sum=$w_peak_sum/$spec_sum;
    }

    if ($den!=0) {
	$w_peak[$k]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
#	$w_peak[$k]=0.5*($w_peak[$k]+$w_peak_sum);
    } else {
	$w_peak[$k]=$w_guess;
    }
$flux_peak[$k]=$spec[$nb];




if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.6);           # Set character height
#    pgenv($jmin,$jmax,$min-5,$max+5,0,0);
    pgenv($Jmin,$Jmax,$min-5,$max+5,0,0);
    pgsch(1.6);           # Set character height
    pglabel("X-axis","Counts","X=$x, $k/$naxis[1]");
    pgsch(2.2);           # Set character height
    pgsci(2);
    pgline($naxis[0],\@wave,\@Spec);    
    pgsci(1);
    pgline($naxis[0],\@wave,\@spec);    
    pgsci(7);
    $x=$i_peak[$k];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],1);
    pgsci(4);
    $x=$w_peak[$k];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],2);
    pgsci(7);
#    pgptxt([$x],[0.5*$y],1,0,"$x");
    pgsci(1);
    pgclose;
    pgend;
#    print "$x $y\n";
#    print "Press Enter"; <stdin>;
	if (($command ne "A")&&($command ne "0")) {
	    print "Press Enter (A=Automatic, 0=No plot)";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }
	    if ($command eq "N") {
		$plot=2;
	    }
	}

}

}

#for ($k=($n_fib+1);$k<$naxis[1];$k++) {
for ($k=($n_fib-1);$k>-1;$k--) {
#    $w_guess=$w_peak[$k+1];
$min=1e12;
$max=-1e12;
    $jmin=$W_guess-$delta_n;
    $jmax=$W_guess+$delta_n;

    for ($nm=0;$nm<$naxis[0];$nm++) {    
	$spec[$nm]=$in_cube[$k][$nm];
	$wave[$nm]=$nm;
	if (($nm>=$jmin)&&($nm<=$jmax)) {
	    if ($max<$spec[$nm]) {
		$max=$spec[$nm];
	    }
	    if ($min>$spec[$nm]) {
		$min=$spec[$nm];
	    }
	}
#	print "$wave[$nm] $median[$nm]\n";
    }

#    $jmin=$w_guess-2*$nsearch;
#    $jmax=$w_guess+2*$nsearch;
    $j=$jmin;
    $peak=0;
    $max_peak=-1e12;
    $i_peak[$k]=$w_guess;
    while ($j<$jmax) {
	if ($nsearch!=0) {
	    $peak=1;
	    for ($i=0;$i<$nsearch;$i++) {
		if ($spec[$j-$i]<$spec[$j-$i-1]) {
		    $peak=0;
		}
		if ($spec[$j+$i]<$spec[$j+$i+1]) {
		    $peak=0;
		}
	    } 
	    if (($peak==1)&&($max_peak<$spec[$j])) {
		$i_peak[$k]=$j;    
		$max_peak=$spec[$j];
	    }
	} else {
	    if ($max_peak<$spec[$j]) {
		$i_peak[$k]=$j;    
		$max_peak=$spec[$j];
	    }
	}
	$j++;
	
    }
$na=$i_peak[$k]-1;
$nb=$i_peak[$k];
$nc=$i_peak[$k]+1;
$a=$wave[$na];
$b=$wave[$nb];
$c=$wave[$nc];
$fa=-$spec[$na];
$fb=-$spec[$nb];
$fc=-$spec[$nc];
$den=($fc-2*$fb+$fa);
if ($den!=0) {
    $w_peak[$k]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
} else {
#    $w_peak[$k]=0;
    $w_peak[$k]=$w_guess;
}
$flux_peak[$k]=$spec[$nb];
#print "$k/$naxis[1] $i_peak[$k] $w_peak[$k]\n";

#print "NM=$nm\n";
if ($plot==1) {
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.6);           # Set character height
#    pgenv($jmin,$jmax,$min-5,$max+5,0,0);
    pgenv($Jmin,$Jmax,$min-5,$max+5,0,0);
    pgsch(1.6);           # Set character height
    pglabel("X-axis","Counts","$k/$naxis[1]");
#    pglabel("X-axis","Counts","Spectral Section");
    pgsch(2.2);           # Set character height
    pgsci(2);
    pgline($naxis[0],\@wave,\@Spec);    
    pgsci(1);
    pgline($naxis[0],\@wave,\@spec);    
    pgsci(7);
    $x=$i_peak[$k];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],1);
    pgsci(4);
    $x=$w_peak[$k];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],2);
    pgsci(1);
    pgclose;
    pgend;
#    print "Press Enter"; <stdin>;
	if (($command ne "A")&&($command ne "0")) {
	    print "Press Enter (A=Automatic, 0=No plot)";
	    $command=<stdin>;
	    chop($command);
	    if ($command eq "Q") {
		$plot=0;
	    }
	    if ($command eq "0") {
		$plot=0;
	    }
	    if ($command eq "N") {
		$plot=2;
	    }
	}

}

}

$min=1e12;
$max=-1e12;
$shift_min=1e12;
$shift_max=-1e12;

$max=-1e12;
open(OUT,">check_eline.all");
print OUT "# SHIFT\n";

my @W_peak;
my @Ashift;
my $NN=0;
for ($j=0;$j<$naxis[1];$j++) {
    $index[$j]=$j;
    $shift[$j]=$w_peak[$j]-$w_guess;
    if ($shift_max<$shift[$j]) {
	$shift_max=$shift[$j];
    }
    if ($shift_min>$shift[$j]) {
	$shift_min=$shift[$j];
    }
    print OUT "$shift[$j]\n";
    $ashift[$j]=abs($shift[$j]);
    if ($flux_peak[$j]>0) {
	$W_peak[$NN]=$w_peak[$j];
	$Ashift[$NN]=$ashift[$j];
	$NN++;
    }
#    print "$w_peak[$j] $flux_peak[$j]\n";
}
close(OUT);

$med_w_peak=median(@W_peak);
$mean_w_peak=mean(@W_peak);
$sig_w_peak=median(@Ashift);

$med_flux_peak=median(@flux_peak);
$mean_flux_peak=mean(@flux_peak);
$sig_flux_peak=sigma(@flux_peak);

if ($plot==1) {
    print "$med_w_peak $mean_w_peak $sig_w_peak $med_flux_peak $mean_flux_peak $sig_flux_peak\n";



}
open(OUT,">check_eline.out");
print OUT "$med_w_peak $mean_w_peak $sig_w_peak $med_flux_peak $mean_flux_peak $sig_flux_peak\n";
close(OUT);

$med_w_peak=median(@w_peak_mean);
$mean_w_peak=mean(@w_peak_mean);

#if ($plot==1) {
#    print "$med_w_peak $mean_w_peak \n";
#
#}

exit;
