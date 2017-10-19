#!/usr/bin/perl

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


if ($#ARGV<2) {
    print "USE:get_index.pl SPECTRUM.TXT REDSHIFT DEV\n";    
    exit;
}

$spec_file=$ARGV[0];
$z=$ARGV[1];
$dev=$ARGV[2];

$min=1e12;
$max=-1e12;
$n=0;
open(FH,"<$spec_file");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $w[$n]=$data[1];
    $f[$n]=$data[2];
    if ($min>$f[$n]) {
	$min=$f[$n];
    }
    if ($max<$f[$n]) {
	$max=$f[$n];
    }
    $n++;
}
close(FH);
#print "$n spectral pixels\n";
$median=median(@f);
$min=-0.1*abs($median);
$max=3.5*$median;


$w_min=1e12;
$w_max=-1e12;
$ni=0;
while($line=<DATA>) {
    chop($line);
    @data=split(" ",$line);
    $name[$ni]=$data[0];
    $L1[$ni]=$data[1]*(1+$z);
    $L2[$ni]=$data[2]*(1+$z);
    $Lb1[$ni]=$data[3]*(1+$z);
    $Lb2[$ni]=$data[4]*(1+$z);
    $Lr1[$ni]=$data[5]*(1+$z);
    $Lr2[$ni]=$data[6]*(1+$z);
    $Lb[$ni]=($Lb1[$ni]+$Lb2[$ni])/2;
    $Lr[$ni]=($Lr1[$ni]+$Lr2[$ni])/2;
    if ($Lb1[$ni]<$w_min) {
	$w_min=$Lb1[$ni];
    }
    if ($Lr2[$ni]>$w_max) {
	$w_max=$Lr2[$ni];
    }
    $ni++;
}
#print "$ni indeces added\n";
$w_min=$w_min-200;
$w_max=$w_max+200;

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.6);           # Set character height
#pgenv($w[0],$w[$n-1],$min-5,$max+5,0,0);
pgenv($w_min,$w_max,$min-0.5,$max+0.5,0,0);
pgsch(1.6);           # Set character height
pglabel("Wavelength","Flux","");
pgsch(2.2);           # Set character height
pgsci(1);
#pgline($n,\@w,\@f);    
pgbin($n,\@w,\@f,1);    

print "$spec_file ";
$cdelt=$w[1]-$w[0];
for ($j=0;$j<$ni;$j++) {
    print "$name[$j] ";
	if ($name[$j] !~ "D4000") {

    $Sb=0;
    $Sr=0;
    $nr=0;
    $nb=0;
    $i1_b=($Lb1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
    $i2_b=($Lb2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
    for ($i=int($i1_b+1);$i<int($i2_b);$i++) {
	$Sb=$Sb+$f[$i]*$cdelt;
	$nb++;
    }
    $i=int($i1_b);
    $ff=$i+1-$i1_b;
    $Sb=$Sb+($f[$i])*$ff*$cdelt;
    $i=int($i2_b);
    $ff=$i2_b-$i;
    $Sb=$Sb+($f[$i])*$ff*$cdelt;
    $Sb=$Sb/($Lb2[$j]-$Lb1[$j]);

    $i1_r=($Lr1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
    $i2_r=($Lr2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
    for ($i=int($i1_r+1);$i<int($i2_r);$i++) {
	$Sr=$Sr+$f[$i]*$cdelt;
	$nr++;
    }
    $i=int($i1_r);
    $ff=$i+1-$i1_r;
    $Sr=$Sr+($f[$i])*$ff*$cdelt;
    $i=int($i2_r);
    $ff=$i2_r-$i;
    $Sr=$Sr+($f[$i])*$ff*$cdelt;
    $Sr=$Sr/($Lr2[$j]-$Lr1[$j]);



    $EW=0;
    $k=0;
    $i1=($L1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
    $i2=($L2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
    for ($i=int($i1+1);$i<int($i2);$i++) {
	$C=$Sb*(($Lr[$j]-$w[$i])/($Lr[$j]-$Lb[$j]))+$Sr*(($w[$i]-$Lb[$j])/($Lr[$j]-$Lb[$j]));
	$EW=$EW+(1-$f[$i]/$C)*($w[$i]-$w[$i-1]);
	$CK[$k]=$C;
	$wk[$k]=$w[$i];
	$k++;
    }
    $i=int($i1);
    $ff=$i+1-$i1;
    $C=$Sb*(($Lr[$j]-$w[$i])/($Lr[$j]-$Lb[$j]))+$Sr*(($w[$i]-$Lb[$j])/($Lr[$j]-$Lb[$j]));
    $EW=$EW+(1-$f[$i]/$C)*($w[$i]-$w[$i-1])*$ff;
    $i=int($i2);
    $ff=$i2-$i;
    $C=$Sb*(($Lr[$j]-$w[$i])/($Lr[$j]-$Lb[$j]))+$Sr*(($w[$i]-$Lb[$j])/($Lr[$j]-$Lb[$j]));
    $EW=$EW+(1-$f[$i]/$C)*($w[$i]-$w[$i-1])*$ff;

    $EW=$EW/(1+$z);

	} else {
 $Sb=0;
		    $nb=0;
		    $i1_b=($Lb1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
		    $i2_b=($Lb2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
		    for ($i=int($i1_b+1);$i<int($i2_b);$i++) {
			$Sb=$Sb+$f[$i]*$cdelt;
			$nb++;
		    }
		    $i=int($i1_b);
		    $ff=$i+1-$i1_b;
		    $Sb=$Sb+($f[$i])*$ff*$cdelt;
		    $i=int($i2_b);
		    $ff=$i2_b-$i;
		    $Sb=$Sb+($f[$i])*$ff*$cdelt;
		    $Sb=$Sb/($Lb2[$j]-$Lb1[$j]);

		    $S=0;
		    $k=0;
		    $i1=($L1[$j]-$w[0]-0.5*$cdelt)/$cdelt;
		    $i2=($L2[$j]-$w[0]-0.5*$cdelt)/$cdelt;
		    for ($i=int($i1+1);$i<int($i2);$i++) {
			$S=$S+$f[$i]*$cdelt;
			$k++;
		    }
		    $i=int($i1);
		    $ff=$i+1-$i1;
		    $S=$S+($f[$i])*$ff*$cdelt;
		    $i=int($i2);
		    $ff=$i2-$i;
		    $S=$S+($f[$i])*$ff*$cdelt;
		    $S=$S/($L2[$j]-$L1[$j]);
		    
		    if ($Sb!=0) {
			$EW=$S/$Sb;
		    } else {
			$EW=1e16;
		    }

	}

#    for ($i=1;$i<$n;$i++) {
#	if (($w[$i]>=$L1[$j])&&($w[$i]<=$L2[$j])) {
#	}
#    }
    print " $EW ";
    pgsls(2);
    pgsci(5);
    pgsfs(2);
   pgpoly(4,[$Lb1[$j],$Lb1[$j],$Lb2[$j],$Lb2[$j]],[$max*0.5,$max*0.8,$max*0.8,$max*0.5]);
#    pgline(2,[$Lb2[$j],$Lb2[$j]],[$max*0.5,$max*0.8]);
    pgsci(2);
   pgpoly(4,[$Lr1[$j],$Lr1[$j],$Lr2[$j],$Lr2[$j]],[$max*0.5,$max*0.8,$max*0.8,$max*0.5]);
#    pgline(2,[$Lr1[$j],$Lr1[$j]],[$max*0.5,$max*0.8]);
#    pgline(2,[$Lr2[$j],$Lr2[$j]],[$max*0.5,$max*0.8]);
    pgsls(1);
    pgsci(8);
    pgline($k,\@wk,\@CK);    
    pgpoint(1,[$Lb[$j]],[$Sb],2);
    pgpoint(1,[$Lr[$j]],[$Sr],2);
    pgsci(3);
   pgpoly(4,[$L1[$j],$L1[$j],$L2[$j],$L2[$j]],[$max*0.5,$max*0.8,$max*0.8,$max*0.5]);
    pgsci(1);
}





pgclose;
pgend;
print "\n";


__DATA__
Hd     4083.500 4122.250 4041.600 4079.750 4128.500 4161.000
Hb     4847.875 4876.625 4827.875 4847.875 4876.625 4891.625
Mgb    5160.125 5192.625 5142.625 5161.375 5191.375 5206.375 
Fe5270 5245.650 5285.650 5233.150 5248.150 5285.650 5318.150
Fe5335 5312.125 5352.125 5304.625 5315.875 5353.375 5363.375
D4000  4050.000 4250.000 3750.000 3950.000 0.000    1.000
Hdmod  4083.500 4122.250 4079     4083     4128.500 4161.000
Hg     4319.75  4363.50  4283.50  4319.75  4367.25  4419.75
