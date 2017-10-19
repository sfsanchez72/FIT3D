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


if ($#ARGV<1) {
    print "USE: plot_Age_Met.pl fit_rss_eline.out DEVICE [MIN MAX]\n";
    exit;
}

$infile=$ARGV[0];
$device=$ARGV[1];


#$infile="spec_A2218_".$id.".fit_template";
#$infile=$ARGV[0];



$n=0;
open(FH,"<$infile");
$min=1e12;
$max=-1e-12;
while($line=<FH>) {
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
close(FH);
print "MINMAX = $min $max\n";
if ($#ARGV==3) {
    $min=$ARGV[2];
    $max=$ARGV[3];
}


@age_sort = sort { $a <=> $b } @age;
#@met_sort = sort { $a <=> $b } @met;

$nx=($#age+1)/6;
$k=0;
for ($j=0;$j<6;$j++) {
    $k=0;
    for ($i=0;$i<=$#age;$i=$i+6) {
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

$ny=6;
#print "$nx $ny\n";
@tr=(-0.5,1,0,-0.5,0,1);
pgbeg(0,$device,1,1);
pgpap(10.,0.55);
pgsch(1.2);
#pgenv(0,$nx,0,6,1,-1);
pgenv(0,$nx,0,6,1,-1);

($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
$nc=$pl->getdim(0);
for ($j=0;$j<$nc;$j++) {
    $l[$j] = $pl->slice($j)->sclr;
    $r[$j] = $pr->slice($j)->sclr;
    $g[$j] = $pg->slice($j)->sclr;
    $b[$j] = $pb->slice($j)->sclr;
}

pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contra);



pgsch(1.4);
pglabel("Age (Gyr)","Metallicity","");
pgimag(\@chisq,$nx,$ny,1,$nx,1,$ny,$max*1.5,$min*0.5,\@tr);
#$n=10;
#for ($i=0;$i<$n;$i++) {
#    $plevel[$i]=10**($level[$i]);
#    
#}
#$le
#pgslw($width);
#pgcons(\@a_cont,$n1,$n2,1,$n1,1,$n2,\@plevel,$n,\@tr);
#
#pgcons(\@chisq,$nx,$ny,1,$nx,1,$ny,\@l,$n,\@tr);
pgsci(1);
pgslw(1);
#pgsch(0.9);
pgsch(1.5);
$k=0.5;
for ($i=0;$i<=$#age;$i=$i+6) {
    pgptxt($k,-0.4,0,0.6,"$age_sort[$i]");
    $k=$k+1;
}
$k=0.5;
for ($i=0;$i<6;$i++) {
    pgptxt(-0.2,$k,90,0.5,"$met[$i]");
    $k=$k+1;
}
#print "$x_min $y_min\n";

pgsci(0);
@chi_list_sort = sort { $a <=> $b } @chi_list;
for ($i=0;$i<30;$i++) {
    $xx=$x{$chi_list_sort[$i]};
    $yy=$y{$chi_list_sort[$i]};
    $min_cut=int($chi_list_sort[$i]*100)/100;
    #pgptxt($xx+0.25,$yy+0.5,0,0,"$min_cut");
    $level[$i]=$chi_list_sort[$i];
}

$age_w=0;
$met_w=0;
$age2_w=0;
$met2_w=0;
$agemet_w=0;
$w=0;
$mo=0;
for ($i=0;$i<$n;$i++) {
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
$age_w=$age_w/$w;
$met_w=$met_w/$w;
$age2_w=$age2_w/$w-($age_w)**2;
$met2_w=$met2_w/$w-($met_w)**2;
$agemet_w=$agemet_w/$w-$age_w*$met_w;
$THETA=0.5*atan(2*$agemet_w/($age2_w-$met2_w));
$A=sqrt(0.5*($age2_w+$met2_w)+sqrt((0.5*($age2_w-$met2_w))**2+$agemet_w**2));
$B=sqrt(0.5*($age2_w+$met2_w)-sqrt((0.5*($age2_w-$met2_w))**2+$agemet_w**2));

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
for ($i=0;$i<=($#age-6);$i=$i+6) {    
    if (($age_w>=$age_sort[$i])&&($age_w<$age_sort[$i+6])) {
	$xx=$kk+($age_w-$age_sort[$i])/($age_sort[$i+6]-$age_sort[$i]);
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
print "AGEMET $age_w $sig_age $met_w $sig_met\n"; #$sig_age_w $sig_met_w\n";
pgsci($color);
#pgcons(\@chisq,$nx,$ny,1,$nx,1,$ny,[$level[8]],1,\@tr);
pgsfs(3);
#pgconf(\@chisq,$nx,$ny,1,$nx,1,$ny,0.0,$level[8],\@tr);
pgsch(3);
pgpoint(1,[$xx],[$yy],$sym);
pgsch(1.5);

pgsci(1);
pgbox("SBCI",0,0,"SBCI",0,0);
pgclos();
pgend();

exit;
