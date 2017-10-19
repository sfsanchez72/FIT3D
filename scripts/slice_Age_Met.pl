#!/usr/bin/perl
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
use PDL::Slatec;
use PDL::Image2D;
use PDL::Graphics::LUT;

#use PDL::Matrix;
$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

$ENV{'PGPLOT_ENVOPT'}="I";

#$table="idl5"; $reverse=1; $bright=0.65; $contra=0.55; $color_cont=0;
#$table="ramp"; $reverse=0; $bright=0.65; $contra=0.55; $color_cont=0;
$table="heat"; $reverse=1; $bright=0.65; $contra=0.55; $color_cont=0;

$f=1/cos((66.2117906/180)*3.14159);


if ($#ARGV<3) {
    print "USE: slice_Age_Met.pl fit_rss_eline.out slice.txt n_met PREFIX\n";
    exit;
}

$infile=$ARGV[0];
$slice=$ARGV[1];
$n_met=$ARGV[2];
$prefix=$ARGV[3];
#$device=$ARGV[4];
#$infile="spec_A2218_".$id.".fit_template";
#$infile=$ARGV[0];

open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	if ($x_min>$x[$ns]) {
	    $x_min=$x[$ns];
	}
	if ($y_min>$y[$ns]) {
	    $y_min=$y[$ns];
	}
	if ($x_max<$x[$ns]) {
	    $x_max=$x[$ns];
	}
	if ($y_max<$y[$ns]) {
	    $y_max=$y[$ns];
	}
	if ($ns>0) {
	    if (($dpix>abs($x[$ns]-$x[$ns-1]))&&(abs($x[$ns]-$x[$ns-1])!=0)) {
		$dpix=abs($x[$ns]-$x[$ns-1]);
	    }
	}
	$ns++;
    } else {
	$header=$line;
    }
}
close(FH);



$n=0;
open(FH,"<$infile");
$min=1e12;
$max=-1e-12;
$n_spec=0;
while($line=<FH>) {
    chop($line);
    if ($line =~ "RSS") {
	@data=split(" ",$line);
	$n_now=$data[1];
	$n_mod=$data[3];
    }
    $n=0;
    my @chi_list;
    my @age;
    my @met;
    for ($mod=0;$mod<$n_mod;$mod++) {
#	print "$mod/$n_mod\n";
	$line=<FH>;
	if ($line =~ "ssp") {
	    chop($line);
	    @data=split(" ",$line);
	    $name=$data[0];
	    $name =~ s/\.spec//;
	    $name =~ s/spec_ssp_//;
	    ($AGE,$MET)=split("_",$name);
	    $AGE_NAME[$n]=$AGE;
	    $MET_NAME[$n]=$MET;
	    if ($AGE =~ "Myr") {
		$age[$n]=$AGE;
		$age[$n] =~ s/Myr//;
		$age[$n]=$age[$n]/1000;
	    } else {
		$age[$n]=$AGE;
		$age[$n] =~ s/Gyr//;
	    }
	    $met[$n]=$MET;
	    $met[$n] =~ s/z/0\./;
	    $val=$age[$n]."_".$met[$n];
	    $chisq_val{$val}=$data[1];
	    if ($min>$data[1]) {
		$min=$data[1];
	    }
	    if ($max<$data[1]) {
		$max=$data[1];
	    }
	    $chi_list[$n]=$data[1];
#	print "$name $age[$n] $met[$n] $chisq[$n]\n";
	    $n++;
	}
    }    
#   
    
    @age_sort = sort { $a <=> $b } @age;
#@met_sort = sort { $a <=> $b } @met;
    
    $nx=($#age+1)/$n_met;
#    print "NX=$nx NAGE=$#age NMET=$n_met\n";

    $k=0;
    for ($j=0;$j<$n_met;$j++) {
	$k=0;
	for ($i=0;$i<=$#age;$i=$i+$n_met) {
	    $val=$age_sort[$i]."_".$met[$j];
	    $chisq[$j][$k]=$chisq_val{$val};
	    $x{$chisq_val{$val}}=$k;
	    $y{$chisq_val{$val}}=$j;
	    if ($min==$chisq[$j][$k]) {
		$x_min=$k;
		$y_min=$j;
	    }

	    #print "$val $k $j $chisq[$j][$k]\n";
	    $k++;
	}
    }
    
    $chisq[0][0]=$chisq[0][1];

    $ny=$n_met;
#print "$nx $ny\n";
    @tr=(-0.5,1,0,-0.5,0,1);
    @chi_list_sort = sort { $a <=> $b } @chi_list;
    my @level;
    for ($i=0;$i<$n_mod;$i++) {
	$xx=$x{$chi_list_sort[$i]};
	$yy=$y{$chi_list_sort[$i]};
	$min_cut=int($chi_list_sort[$i]*100)/100;
	#pgptxt($xx+0.25,$yy+0.5,0,0,"$min_cut");
	$level[$i]=$chi_list_sort[$i];
#	print "LEV=$i $level[$i]\n";
    }
    
    $age_w=0;
    $met_w=0;
    $age2_w=0;
    $met2_w=0;
    $agemet_w=0;
    $w=0;
    $mo=0;
    for ($i=0;$i<$n;$i++) {
#	print "TEST='$i' '$chi_list[$i]' '$level[5]'\n";
	if ($chi_list[$i]<$level[5]) {
	    $age_w=$age_w+$age[$i]/$chi_list[$i];
	    $met_w=$met_w+$met[$i]/$chi_list[$i];
	    $age2_w=$age2_w+($age[$i]*$age[$i])/$chi_list[$i];
	    $met2_w=$met2_w+($met[$i]*$met[$i])/$chi_list[$i];
	    $agemet_w=$agemet_w+($age[$i]*$met[$i])/$chi_list[$i];
	    $w=$w+1/$chi_list[$i];
	    $age_mo[$mo]=$age[$i];
	    $met_mo[$mo]=$met[$i];
	    $age_mo_w[$mo]=$age[$i]/$chi_list[$i];
	    $met_mo_w[$mo]=$met[$i]/$chi_list[$i];
	    $mo++;
	}
    }
#	print "W=$w\n";
    if ($w==0) {
	$w=1;
    }
    $age_w=$age_w/$w;
    $met_w=$met_w/$w;
    $age2_w=$age2_w/$w-($age_w)**2;
    $met2_w=$met2_w/$w-($met_w)**2;
    $agemet_w=$agemet_w/$w-$age_w*$met_w;
#    $THETA=0.5*atan(2*$agemet_w/($age2_w-$met2_w));
#    $A=sqrt(abs(0.5*($age2_w+$met2_w)+sqrt((0.5*($age2_w-$met2_w))**2+$agemet_w**2)));
#    $B=sqrt(abs(0.5*($age2_w+$met2_w)-sqrt((0.5*($age2_w-$met2_w))**2+$agemet_w**2)));
    
    $var_age=0;
    $var_met=0;
    $var_age_w=0;
    $var_met_w=0;
    for ($i=0;$i<$mo;$i++) {
	$var_age=$var_age+($age_w-$age_mo[$mo])**2;
	$var_met=$var_met+($met_w-$met_mo[$mo])**2;
	$var_age_w=$var_age_w+($age_w*$w-$age_mo_w[$mo])**2;
	$var_met_w=$var_met_w+($met_w*$w-$met_mo_w[$mo])**2;
#    print "$age_w $age_mo[$i]\n";
    }
    $sig_age=0.5*sqrt($var_age)/($mo-1);
    $sig_met=sqrt($var_met)/($mo-1);
    $sig_age_w=(0.5*sqrt($var_age_w)/($mo-1))/$w;
    $sig_met_w=(sqrt($var_met_w)/($mo-1))/$w;
    
    $xx=0;
    $kk=0;
    for ($i=0;$i<=($#age-$n_met);$i=$i+$n_met) {    
	if (($age_w>=$age_sort[$i])&&($age_w<$age_sort[$i+$n_met])) {
	    $xx=$kk+($age_w-$age_sort[$i])/($age_sort[$i+$n_met]-$age_sort[$i]);
	}
	$kk++;
    }
    for ($i=0;$i<5;$i++) {
	if (($met_w>=$met[$i])&&($met_w<$met[$i+1])) {
	    $yy=$i+($met_w-$met[$i])/($met[$i+1]-$met[$i]);
	}
    }
    $xx=$xx+0.25;
    $yy=$yy+0.25;
#$xx=$x{$chi_list_sort[0]}+0.25;
#$yy=$y{$chi_list_sort[0]}+0.5;
    if ($nser>2.5) {
	$color=2;
	$sym=17;
    } else {
	if ($nser>1.75) {
	    $color=8;
	    $sym=16;
	} else {
	    $color=3;
	    $sym=3;
	}
    }
    

#print "$id $age_w $met_w $A $B $THETA $sig_age $sig_met\n";
#    print "AGEMET $age_w $sig_age $met_w $sig_met\n"; #$sig_age_w $sig_met_w\n";

	$AGE_W[$n_spec]=$age_w;
	$S_AGE_W[$n_spec]=$sig_age;
	$MET_W[$n_spec]=$met_w;
	$S_MET_W[$n_spec]=$sig_met;
    
    $n_spec++;
}
close(FH);

$age_file=$prefix."_age.txt";
open(FH,">$age_file");
#    print FH "$header\n";
for ($i=0;$i<$ns;$i++) {
    print FH "$id[$i] $x[$i] $y[$i] $AGE_W[$i]\n";
}
close(FH);

$e_age_file=$prefix."_e_age.txt";
open(FH,">$e_age_file");
#print FH "$header\n";
for ($i=0;$i<$ns;$i++) {
    print FH "$id[$i] $x[$i] $y[$i] $S_AGE_W[$i]\n";
}
close(FH);

$met_file=$prefix."_met.txt";
open(FH,">$met_file");
#print FH "$header\n";
for ($i=0;$i<$ns;$i++) {
    print FH "$id[$i] $x[$i] $y[$i] $MET_W[$i]\n";
}
close(FH);

$e_met_file=$prefix."_e_met.txt";
open(FH,">$e_met_file");
#print FH "$header\n";
for ($i=0;$i<$ns;$i++) {
    print FH "$id[$i] $x[$i] $y[$i] $S_MET_W[$i]\n";
}
close(FH);

$pdl_age=zeroes(40,$ns);
for ($j=0;$j<$ns;$j++) {
    for ($i=0;$i<10;$i++) {
	set($pdl_age,$i,$j,$AGE_W[$j]);
    }
    for ($i=10;$i<20;$i++) {
	set($pdl_age,$i,$j,$S_AGE_W[$j]);
    }
    for ($i=20;$i<30;$i++) {
	set($pdl_age,$i,$j,$MET_W[$j]);
    }
    for ($i=30;$i<40;$i++) {
	set($pdl_age,$i,$j,$S_MET_W[$j]);
    }
}
$pdl_file=$prefix.".age_met.fits";
#$pdl_age->
$h=$pdl_age->gethdr;
$$h{CRPIX1}=1;
$$h{CRVAL1}=1;
$$h{CDELT1}=1;
$pdl_age->sethdr($h);
$pdl_age->wfits($pdl_file);

print "Results written in $prefix\* files\n";

exit;
