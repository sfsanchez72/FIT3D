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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<1) {
    print "USE: spec_plot.pl SET_INPUTFILES DEVICE [MIN MAX] [WMIN WMAX] [NAME] [WIDTH]\n";
    exit;
}

$width=5;
$sinput=$ARGV[0];
$dev=$ARGV[1];
$def=0;
if ($#ARGV==3) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $def=1;
}


$y_min=1e12;
$y_max=-1e12;
$nf=0;
open(C,"<$sinput");
while ($line=<C>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==0) {
	$file[$nf]=$data[0];
	$ratio[$nf]=1;
    } else {
	$file[$nf]=$data[0];
	$ratio[$nf]=$data[1];
    }
#    print "$ratio[$nf]\n";
    $nf++;
}
close(C);
#print "$nf files found\n";
$wmin=1e12;
$wmax=-1e12;
for ($i=0;$i<$nf;$i++) {
    $input=$file[$i];
#    print "$input ";
    $n=0;
    open(FH,"<$input");
    while ($line=<FH>) {
	chop($line);
	if (($line !~ "#")&&($line !~ "wave")) {
	    @data=split(" ",$line);
	    if ($#data==2) {
		$flux[$i][$n]=$data[2]*$ratio[$i];	
		$wave[$i][$n]=$data[1];
		if ($wmin>$data[1]) {
		    $wmin=$data[1];
	    }
		if ($wmax<$data[1]) {
		    $wmax=$data[1];
		}
	    } else {
		$flux[$i][$n]=$data[1]*$ratio[$i];	
		$wave[$i][$n]=$data[0];
		
		if ($wmin>$data[0]) {
		    $wmin=$data[0];
		}
		if ($wmax<$data[0]) {
		    $wmax=$data[0];
		}
	    }
#	print "$flux[$i][$n] $wave[$i][$n] $ratio[$i]\n";

	if ($flux[$i][$n]>$y_max) {
	    $y_max=$flux[$i][$n];
	}
	if ($flux[$i][$n]<$y_min) {
	    $y_min=$flux[$i][$n];
	}
	$n++;
	}
    }
    $nline[$i]=$n;
    close(FH);
#    print "$n\n";
}
#print "$n\n";
#$wmin=$wave[0][0];
#$wmax=$wave[0][$n-1];
#print "[$wmin,$wmax]\n";
if ($#ARGV==5) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $def=1;
}

$name="";
if ($#ARGV==6) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $name=$ARGV[6];
    $def=1;
}

if ($#ARGV==7) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $wmin=$ARGV[4];
    $wmax=$ARGV[5];
    $name=$ARGV[6];
    $width=$ARGV[7];
    $def=1;
}




if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}


#@stats=stats(pdl(@flux));
#$y_max=1.5*$stats[2];

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength (\\A)","Flux (10\\u-16\\d Erg s\\u-1\\d cm\\u-2\\d \\A\\u-1\\d )","");
pgsch(1.65);
pgptxt($wmin+0.85*($wmax-$wmin),$y_min+0.9*($y_max-$y_min),0,0,$name);
pgsch(1.4);
$c=1;
for ($i=0;$i<$nf;$i++) {    
    $n=$nline[$i];#=$n;
#    print "$file[$i] $n $wave[$i][0] $wave[$i][$n-1] ";
    $wave_mean=0;
    $sum_flux=0;
    my @wave_now;
    my @flux_now;
    for ($j=0;$j<$n;$j++) {
	$wave_now[$j]=$wave[$i][$j];
	$flux_now[$j]=$flux[$i][$j];
	if ($flux_now[$j]>0) {
	    $wave_mean=$wave_mean+$flux_now[$j]*$wave_now[$j];
	    $sum_flux=$sum_flux+$flux_now[$j];
	}
    }
    if ($sum_flux>0) {
	$wave_mean=$wave_mean/$sum_flux;
#	print "$wave_mean\n";
    }

   
    


    pgsci($c);
    $style=$i+1;
    pgsls($style);
    pgslw(6/($i+1));
    pgline($n,\@wave_now,\@flux_now);
    pgslw(1);
    pgsls(1);
    pgsch(1.8);
    $kk=$i+1;    
    srandom;
    $pos=int($n/5+$i*5);#+random(3);
#    pgsci(1);
    #pgptxt($wave_now[$pos],$flux_now[$pos],0,0,$kk);


    if ($i>0) {
	my @flux_new;
	my @rat_new;

	($wave_min,$wave_max)=minmax(@wave_now);
	$j=0;
	while ($wave_old[$j]<$wave_max) {
	    $j++;
	}
	$n_old=$j-1;

	$out_spec_pdl = interpol(pdl(@wave_old), pdl(@wave_now), pdl(@flux_now));



	my @flux_new_tmp=list($out_spec_pdl);
#	@flux_new=@flux_new;
#	$width=5;
	for ($j=0;$j<$n_old;$j++) {
	    $flux_new[$j]=$flux_old[$j];
	    
	}
	for ($j=2*$width;$j<$n_old-2*$width;$j++) {
	    $flux_new[$j]=$flux_old[$j];
	    $sum=0;
	    $C=0;
	    for ($k=$j-$width;$k<$j+$width;$k++) {
		$sum=$sum+$flux_new_tmp[$k]*exp(-0.5*(($k-$j)/(0.5*$width))**2);
		$C=$C+exp(-0.5*(($k-$j)/(0.5*$width))**2);
	    }
	    
	    $flux_new[$j]=$sum/$C;#-$flux_old[$j];	    
	    
	}

	for ($j=0;$j<$n_old;$j++) {
	    $flux_new[$j]=$flux_new[$j]-$flux_old[$j];
	    $rat_new[$j]=$flux_new[$j]/$flux_old[$j];
	}
	
	pgsci($c+1);
	pgline($n_old,\@wave_old,\@flux_new);

	@stats=stats(pdl(@rat_new));
	print "$name @stats\n";

	open(OUT,">comp_spectra_res.out");
	for ($i=0;$i<$n_old;$i++) {
	    if (($wave_old[$i]>$wmin)&&($wave_old[$i]<$wmax)) {
		if (($rat_new[$i]>-1)&&($rat_new[$i]<1)) {
		    print OUT "$i $wave_old[$i] $rat_new[$i]\n";
		}
	    }
	}
	close(OUT);
	

    }
    pgsls(1);
    pgsch(1.2);



    $c++;
    if ($c>15) {
	$c=1;
    }
    @wave_old=@wave_now;
    @flux_old=@flux_now;
    $n_old=$n;

}
pgclose;
pgend;



exit;
