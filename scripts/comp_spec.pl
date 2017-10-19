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


if ($#ARGV<2) {
    print "USE: comp_spec.pl SPEC1.txt SPEC2.txt ratio.txt [wave_min wave_max]\n [NSMOOTH]";
    exit;
}

$file1=$ARGV[0];
$file2=$ARGV[1];
$out_file=$ARGV[2];
$wave_min=-1e12;
$wave_max=1e12;
if ($#ARGV==4) {
    $wave_min=$ARGV[3];
    $wave_min=$ARGV[4];
}
$NS=0;
if ($#ARGV==5) {
    $wave_min=$ARGV[3];
    $wave_min=$ARGV[4];
    $NS=$ARGV[5];
}
$min=1e12;
$max=-1e12;
$min_rat=1e12;
$max_rat=-1e12;
$n1=0;
open(FH,"<$file1");
while($line=<FH>) {
    @data=split(" ",$line);
    $index1[$n1]=$data[0];
    $wave1[$n1]=$data[1];
    $flux1_ini[$n1]=$data[2];
    if ($min>$flux1_ini[$n1]) {
	$min=$flux1_ini[$n1];
    }
    if ($max<$flux1_ini[$n1]) {
	$max=$flux1_ini[$n1];
    }
    $n1++;
}
close(FH);

$n2=0;
$k=0;
open(FH,"<$file2");
while($line=<FH>) {
    @data=split(" ",$line);
    $index2[$n2]=$data[0];
    $wave2[$n2]=$data[1];
    $flux2_ini[$n2]=$data[2];
    if ($flux1_ini[$n2]>0) {
	$RATIO[$n2]=$flux2_ini[$n2]/$flux1_ini[$n2];
    }

    if (($wave2[$n2]>$wave_min)&&($wave2[$n2]<$wave_max)) {
	$ratio[$k]=$flux2_ini[$n2]/$flux1_ini[$n2];
	$diff[$k]=$flux2_ini[$n2]-$flux1_ini[$n2];
	$k++;
    }

    if ($min>$flux2_ini[$n2]) {
	$min=$flux2_ini[$n2];
    }
    if ($max<$flux2_ini[$n2]) {
	$max=$flux2_ini[$n2];
    }
    if ($min_rat>$RATIO[$n2]) {
	$min_rat=$RATIO[$n2];
    }
    if ($max_rat<$RATIO[$n2]) {
	$max_rat=$RATIO[$n2];
    }
    $n2++;
}
close(FH);
open(OUT,">comp_spec.out");
$median=median(@ratio);
$sigma=sigma_m(@ratio);
print "ratio=$median+-$sigma\n";
print OUT "$median $sigma\n";
$median=median(@diff);
$sigma=sigma_m(@diff);
print "diff=$median+-$sigma\n";
print OUT "$median $sigma\n";
close(OUT);

pgbegin(0,"666/xs",1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
pgsubp(1,2);
pgenv($wave1[0],$wave1[$n1-1],$min,$max,0,0);
pgsch(1.2);           # Set character height
pglabel("Wavelength","Counts","");
pgsci(1);
pgpoint($n1,\@wave1,\@flux1_ini,5);    
pgsci(2);
pgpoint($n2,\@wave2,\@flux2_ini,5);    
pgsci(1);
pgenv($wave1[0],$wave1[$n1-1],$min_rat,$max_rat,0,0);
pglabel("Wavelength","Ratio","");
pgsch(0.9);           # Set character height
pgsci(3);
pgpoint($n1,\@wave1,\@RATIO,2);       
pgsci(2);
#$w1=$wave1[0];
#$w2=$wave1[$n1-1];
pgline(2,[$wave1[0],$wave1[$n1-1]],[$median,$median]);
pgsci(8);
pgsls(2);
pgline(2,[$wave1[0],$wave1[$n1-1]],[$median-$sigma,$median-$sigma]);
pgline(2,[$wave1[0],$wave1[$n1-1]],[$median+$sigma,$median+$sigma]);
pgsls(1);
pgsci(1);
pgclose;
pgend;


if ($NS>0) {
    @smooth=median_filter($NS,\@RATIO);
}


open(FH,">$out_file");
for ($j=0;$j<$n1;$j++) {
    if ($NS>0) {
	$RATIO[$j]=$smooth[$j];
    }
    print FH "$j $wave1[$j] $RATIO[$j]\n";
}
close(FH);

exit;

