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
use PDL::Core;
use PDL::Graphics::LUT;
use Carp;



$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<13) {
    print "USE: maps_orion.pl map.fits min max bright contrast Label factor dev SYSTEMIC_VEL MAP_CONT NLEVELS MIN STEP NWHITE[49] [ROT]\n";
    exit;
}
$mapfile=$ARGV[0];
$min=$ARGV[1];
$max=$ARGV[2];
$bright=$ARGV[3];
$contrast=$ARGV[4];
$label=$ARGV[5];
$factor=$ARGV[6];
$dev=$ARGV[7];
$vw=$ARGV[8];
$map_cont=$ARGV[9];
$nlevels=$ARGV[10];
$min_cont=$ARGV[11];
$steps=$ARGV[12];
$NW=$ARGV[13];

if ($dev !~ "TPNG") { 
    $ENV{PGPLOT_FOREGROUND} = "black";
    $ENV{PGPLOT_BACKGROUND} = "white";
} else {
    $ENV{PGPLOT_FOREGROUND} = "white";
    $ENV{PGPLOT_BACKGROUND} = "white";
}



$rot=0;
if ($#ARGV==14) {
     $rot=$ARGV[14];
}

$color_cont=0;

#$pdl_mask=rfits("M74_mos.mask.fits");
#$pdl_mask=rfits("mask_sc.fits");
#$pdl_map=rfits($mapfile);
#$pdl_cont=rfits($map_cont);



$pdl_map_tmp=rfits($mapfile);
$pdl_cont_tmp=rfits($map_cont);

$is_zero=sum(abs($pdl_map_tmp));



if ($rot==0) {
    $pdl_map=$pdl_map_tmp;
    $pdl_cont=$pdl_cont_tmp;
} else {
#    $pdl_map=$pdl_map_tmp->rot2d($rot,0,1);
#    $pdl_cont=$pdl_cont_tmp->rot2d($rot,0,1);
    $pdl_map=transpose($pdl_map_tmp);
    $pdl_cont=transpose($pdl_cont_tmp);
}



#$crval1=$pdl_map->hdr->{CRVAL1};
#$cdelt1=$pdl_map->hdr->{CDELT1};
#$crval2=$pdl_map->hdr->{CRVAL2};
#$cdelt2=$pdl_map->hdr->{CDELT2};
#print "$crval1 $cdelt1 $crval2 $cdelt2\n";
($nx,$ny)=$pdl_map->dims;
$h=$pdl_map->gethdr;
$crval1=-138.839996337891;
$cdelt1=1;
$crval2=-210.726806640625;
$cdelt2=1;


$pdl_map=$pdl_map*$factor;


$NNx=$nx-1;
$NNy=$ny-1;
$rpdl_map=$pdl_map->slice("$NNx:0,$NNy:0");
@map=list($rpdl_map);
$rpdl_cont=$pdl_cont->slice("$NNx:0,$NNy:0");
@cont=list($rpdl_cont);
$crval1=-30;
$crval2=-40;
$x_max=$crval1;
$x_min=$crval1+$cdelt1*$nx;
$y_min=$crval2;
$y_max=$crval2+$cdelt2*$ny;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	if ($map[$i+$j*$nx]>0) {
	    $map[$i+$j*$nx]=$map[$i+$j*$nx]-$vw;
	}

	if ($map[$i+$j*$nx]>1e12) {
	    $map[$i+$j*$nx]=0;
	}
	if ($map[$i+$j*$nx] eq "nan") {
	    $map[$i+$j*$nx]=0;
	}


	if ($map[$i+$j*$nx]!=(1*$map[$i+$j*$nx])) {
	    $map[$i+$j*$nx]=0;
	}
	if ($map[$i+$j*$nx]==$vw) {
	    $map[$i+$j*$nx]=-1e12;
	}
#	print "$map[$i+$j*$nx]\n";

    }
}


@tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
#@tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
#@tr=($x_max,-$cdelt1,0,$crval2,0,$cdelt2);



#@tr=(0,1,0,0,0,1);
#$x_min=0;
#$x_max=$nx;
#$y_min=0;
#$y_max=$ny;

pgbeg(0,$dev,1,1);
pgscf(2.0);
pgsch(1.2);


$nc=256;
$c=0.9;
$b=1/(128+$c*128**2);
for ($j=0;$j<=128;$j++) {
    $l[$j]=$j*(1/256);
    $r[$j]=$b*($j+$c*$j**2);
    $g[$j]=$j*(2/256);
    $b[$j]=1.0;
#    print  "$j $g[$jj]\n";
}
for ($j=128;$j<256;$j++) {
    $l[$j] = $j*(1/256);
    $r[$j] = 1;
    $g[$j] = ($g[256-$j]);
    $b[$j] = ($r[256-$j]);
}

open(CAL,">CALIFA.vel.ctab\n");
print CAL "# 1 index\n";
print CAL "# 2 Red\n";
print CAL "# 3 Green\n";
print CAL "# 4 Blue\n";
print CAL "# 5 Luminisity\n";
for ($j=0;$j<256;$j++) {
    if (abs($j-128)>1) {
	$F=0.85/exp(-(($j-128)/512)**4);    
    } else {
	$F=1;       
    }
    $rr[$j]=$r[$j]*$F;
#    $gg[$j]=$g[$j]*$F;
    $gg[$j]=$g[$j];
    $bb[$j]=$b[$j]*$F;
    print CAL "$j $rr[$j] $gg[$j] $bb[$j] $l[$j]\n";
#
}
close(CAL);



#pgctab(\@l,\@r,\@g,\@b,$nc,$c,0.5);

pgctab(\@l,\@rr,\@gg,\@bb,$nc,$c,0.5);

#pgenv($x_min,$x_max,$y_min,$y_max,1,0);
pgenv($x_min,$x_max,$y_min,$y_max,1,0);
#pgpap(6.0,1.0);

#pgsch(0.7);           # Set character height
#pgsch(1.2);

if ($dev !~ "TPNG") { 
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
}
if ($is_zero!=0) {
    pgimag(\@map,$nx,$ny,1,$nx,1,$ny,$min,$max,\@tr);
}
pgsci(0);
pgsch(1);
#pgsch(1.1);
@tmp;
$nt=0;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$map[$i+$j*$nx];
#	if ($map[$i+$j*$nx]<=0) {
	if ($val<=0) {
	    $X=$crval1+$cdelt1*$i+1;
	    $Y=$crval2+$cdelt2*$j+1;
	#    print "$i,$j $X,$Y MAP=$map[$i+$j*$nx]\n";
#	    pgpoint(1,[$X],[$Y],16);
	} else {
	    if ($val>6) {
		$tmp[$nt]=$val;
#		print "$tmp[$nt] $nt\n";
		$nt++;
	    }
	}
    }
}
$mean=mean(@tmp);
$sigma=sigma(@tmp);
$sum=$nt*$mean;
print "FLUX=$sum\n";
print "MEAN=$mean+-$sigma\n";


pgsch(3);
pgsci(3);
pgslw(3);
#pgpoint1,[3.65],[-3.13],3);
#pgpoint(1,0,0,3);
pgsch(1.2);
pgslw(1);


pgsci(1);
pgslw(3);
#$nlevels=10;
#$min_cont=$ARGV[11];
#$steps=$ARGV[12];
for ($i=0;$i<$nlevels;$i++) {    
    $levels[$i]=$min_cont+$steps*($i**1.5);
#    print "$i $levels[$i]\n";
}

pgsci(1);
pgslw(4);

if ($dev !~ "TPNG") { 
    pgsci(1);
    pgslw(4);
} else { 
    pgsci(14);
    pgslw(8);
}

if ($is_zero!=0) {
    pgcont(\@cont,$nx,$ny,1,$nx,1,$ny,\@levels,$nlevels,\@tr);
}
pgslw(1);
pgsci(1);

if ($dev !~ "TPNG") { 
    pgsci(1);
    pgsch(1.4);
    pgptxt(-15,28,0,0.5,$label);
    pgsci(1);


pgsch(1.1);
#pgptxt(-85,160,0,0,$label);
#pgptxt(-110,160,0,0,$label);

for ($j=0;$j<256;$j++) {
#    if (abs($j-128)>2) {
	$F=0.85/exp(-(($j-128)/512)**4);    
#    } else {
#	$F=1;       
#    }
    $rr[$j]=$r[$j]*$F;
    $gg[$j]=$g[$j]*$F;
    $bb[$j]=$b[$j]*$F;
    
}
pgctab(\@l,\@rr,\@gg,\@bb,$nc,$c,0.5);
pgwedg("RIV",0.1,4.5,$min,$max,"");
#pgwedg("RI",0,1.0,1,500,"");
pgbox("SBC",0,0,"SBC",0,0);


}




pgclos();
pgend();



exit;
