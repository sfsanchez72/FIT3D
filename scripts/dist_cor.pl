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
    print "USE: dist_cor.pl EXTRACTED.fits CORRECTED.fits DISTORSION_CORRECTION.txt SMOOTH[0/1] start_index delta_index nsearch plot [n_fib] [nsigma] [center_index]\n";
    exit;
}

$input_cube=$ARGV[0];
$output_cube=$ARGV[1];
$output_ps=$ARGV[1];
$dist_cor_file=$ARGV[2];
$output_ps =~ s/.fits/.ps/g;
$smooth=$ARGV[3];
$w_guess=$ARGV[4];
$delta_n=$ARGV[5];
$nsearch=$ARGV[6];
$plot=$ARGV[7];
$n_fib=0;
$nsigma=2;
if ($#ARGV==8) {
    $n_fib=$ARGV[8];
}
if ($#ARGV==9) {
    $n_fib=$ARGV[8];
    $nsigma=$ARGV[9];
}
$center_index=0;
if ($#ARGV==10) {
    $n_fib=$ARGV[8];
    $nsigma=$ARGV[9];
    $center_index=$ARGV[10];
}
$max_delta=1e12;
if ($#ARGV==11) {
    $n_fib=$ARGV[8];
    $nsigma=$ARGV[9];
    $center_index=$ARGV[10];
    $max_delta=$ARGV[11];
}
#print "$#ARGV C.I.=$center_index\n";



print "Paso $input_cube\n";
$nax=read_naxes($input_cube);   
print "Paso\n";
@naxis=@$nax;
print "$naxis[0] $naxis[1]\n"; 
print "Reading file $input_cube...\n";
#exit;
@in_cube=read_img($input_cube);

print "Readed...\n";
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
print "We start look for the 1st peak ...\n\n";


$j=$jmin;
$peak=0;
$max_peak=-1e12;
#while (($peak==0)&&($j<$jmax)) {
$i_peak[$n_fib]=$w_guess;
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
if ($den!=0) {
    $w_peak[$n_fib]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
} else {
    $w_peak[$n_fib]=$w_guess;
}
if ($center_index==0) {
    $center_index=$w_peak[$n_fib];
}
#print "$i_peak[$n_fib] $w_peak[$n_fib]\n";

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
    $x=$i_peak[$n_fib];
    $y=$spec[$x];
    pgpoint(1,[$x],[$y],1);
    pgsci(4);
    $x=$w_peak[$n_fib];
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

print "We look for the peaks in all the spectra\n";
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
#	print "$wave[$nm] $median[$nm]\n";
    }

#    $jmin=$w_guess-2*$nsearch;
#    $jmax=$w_guess+2*$nsearch;
    $peak=0;
    $max_peak=-1e12;
    $j=$jmin;
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
#$min=1e12;
#$max=-1e12;
$shift_min=1e12;
$shift_max=-1e12;

$max=-1e12;
for ($j=0;$j<$naxis[1];$j++) {
    $index[$j]=$j;
#    $shift[$j]=$w_peak[$j]-$w_peak[$n_fib];
    $shift[$j]=$w_peak[$j]-$center_index;
    if ($shift_max<$shift[$j]) {
	$shift_max=$shift[$j];
    }
    if ($shift_min>$shift[$j]) {
	$shift_min=$shift[$j];
    }

}

#exit;



#
#cleaning
#
my @new_shift=@shift;
my @new_index=@index;
$ks=0;
$nbox=5;
#for ($j=$nbox;$j<$naxis[1]-$nbox;$j++) {
for ($j=0;$j<$naxis[1];$j++) {
    my @tmp;
    $k=0;
    $border=0;
    if ($j<$nbox) {
	$j_min=$j;
	$j_max=$j+2*$nbox;
	$border=1;
    } else {
	if ($j>=($naxis[1]-$nbox)) {
	    $j_min=$j-2*$nbox;
	    $j_max=$j;
	    $border=1;
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
    if (abs($shift[$j]-$me)<$nsigma*$sig) {
	$new_shift[$ks]=$shift[$j];
	$new_index[$ks]=$index[$j];
	$ks++;
    } else {
	if ($border==0) {
	    $shift[$j]=$me;
	}
    }
}




print "DONE\n";
#
# Check the offsets
#
if ($smooth==1) {
    $npoly=3;
    $shift_pdl = pdl(@new_shift);
    $index_pdl = pdl(@new_index);
    ($s_y,$coeff) = fitpoly1d $index_pdl,$shift_pdl,$npoly;
    my @shift_all;
    for ($j=0;$j<$npoly;$j++) {
	$c[$j]=$coeff->slice($j)->sclr;
#    print "$c[$j] ";
}
    for ($j=0;$j<$naxis[1];$j++) {
	$shift_all[$j]=0;
	$index_all[$j]=$j;
	for ($k=0;$k<$npoly;$k++) {
	    $shift_all[$j]=$shift_all[$j]+$c[$k]*($j**$k);
	}
#    print "$j $shift_all[$j]\n";
	
    }
} else {
    for ($j=0;$j<$naxis[1];$j++) {
	$index_all[$j]=$j;
	$shift_all[$j]=$shift[$j];
    }
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
    pgpoint($naxis[1],\@index,\@shift,1);
    pgsci(3);
    pgpoint($naxis[1],\@new_index,\@new_shift,22);
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
    $min=1e12;
    $max=-1e12;
#	print "$j $i $shift_all[$j][$i]\n";
    if ($shift_all[$j]!=0) {
	for ($k=0;$k<$naxis[0];$k++) {    
	    $in_spec[$k]=$in_cube[$j][$k];
	    $in_wave[$k]=$k;
#	    $out_wave[$k]=$k-$shift_all[$j];
	    $out_wave[$k]=$k-$shift_all[$j];
	    if (($k>=$jmin)&&($k<=$jmax)) {
		if ($in_spec[$k]<$min) {
		    $min=$in_spec[$k];
		}
		if ($in_spec[$k]>$max) {
		    $max=$in_spec[$k];
		}
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
#    pgenv($w_guess,$w_guess+$delta_w*$nm,-100,1200,0,0);
#	    pgenv($w_guess,$w_guess+$delta_w*$nm,$min-1,$max+1,0,0);
#	    $min=0.5*$max;

	    if ($max<0) {
		$max=100;
		$min=-100;
	    }
	    pgenv($Jmin,$Jmax,$min-1,$max+500,0,0);
#	    pgenv($start_w,$start_w+$delta_w*$naxis[2],-100,1200,0,0);
	    pgsch(1.6);           # Set character height
	    pglabel("Id.","Counts","($i,$j)");
	    pgsch(2.2);           # Set character height
	    pgsci(1);
	    pgline($nm,\@wave,\@Spec);    
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











