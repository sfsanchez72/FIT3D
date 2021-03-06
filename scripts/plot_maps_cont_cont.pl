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
use PDL::Image2D;

#$ENV{PGPLOT_FOREGROUND} = "black";
#$ENV{PGPLOT_BACKGROUND} = "white";

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<17) {
    print "USE: maps_orion.pl map.fits min max bright contrast Label factor dev VAL_TO_WHITE MAP_CONT NLEVELS MIN STEP  MAP_CONT NLEVELS MIN STEP NWHITE[49] [TABLE REVERSE] [ROT]\n";
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
$map_cont2=$ARGV[13];
$nlevels2=$ARGV[14];
$min_cont2=$ARGV[15];
$steps2=$ARGV[16];
$NW=$ARGV[17];

if ($dev !~ "TPNG") {
    $ENV{PGPLOT_FOREGROUND} = "black";
    $ENV{PGPLOT_BACKGROUND} = "white";
} else {
    $ENV{PGPLOT_FOREGROUND} = "white";
    $ENV{PGPLOT_BACKGROUND} = "white";
}


$table="idl5"; $reverse=1; 
if ($#ARGV==19){
    $table=$ARGV[18];
    $reverse=$ARGV[19];
}
$rot=0;
if ($#ARGV==20){
    $table=$ARGV[18];
    $reverse=$ARGV[19];
    $rot=$ARGV[20];
}




#$table="smooth2"; $reverse=0; 
#$table="real"; $reverse=0; 
#$bright=0.7; 
#$contra=0.5; 
$color_cont=0;

#$pdl_mask=rfits("M74_mos.mask.fits");
#$pdl_mask=rfits("mask_sc.fits");
$pdl_map_tmp=rfits($mapfile);
$pdl_cont_tmp=rfits($map_cont);
$pdl_cont2_tmp=rfits($map_cont2);

$is_zero=sum(abs($pdl_map_tmp));

if ($rot==0) {
    $pdl_map=$pdl_map_tmp;
    $pdl_cont=$pdl_cont_tmp;
    $pdl_cont2=$pdl_cont2_tmp;
} else {
    $pdl_map=transpose($pdl_map_tmp);
    $pdl_cont=transpose($pdl_cont_tmp);
    $pdl_cont2=transpose($pdl_cont2_tmp);
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
$rpdl_cont2=$pdl_cont2->slice("$NNx:0,$NNy:0");
@cont2=list($rpdl_cont2);
$crval1=-30;
$crval2=-40;
$x_max=$crval1;
$x_min=$crval1+$cdelt1*$nx;
$y_min=$crval2;
$y_max=$crval2+$cdelt2*$ny;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
#	if ($map[$i+$j*$nx]<0) {
#	    $map[$i+$j*$nx]=0;
#	}
	if ($map[$i+$j*$nx]>1e12) {
	    $map[$i+$j*$nx]=0;
	}
	if ($map[$i+$j*$nx] eq "nan") {
	    $map[$i+$j*$nx]=0;
#
	}
	if ($map[$i+$j*$nx]!=(1*$map[$i+$j*$nx])) {
#	    print "$i $j $map[$i+$j*$nx]\n";
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
($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
$nc=$pl->getdim(0);
for ($j=0;$j<$nc;$j++) {
    $l[$j] = $pl->slice($j)->sclr;
    $r[$j] = $pr->slice($j)->sclr;
    $g[$j] = $pg->slice($j)->sclr;
    $b[$j] = $pb->slice($j)->sclr;
#    print "$j/$nc $l[$j] $r[$j] $g[$j] $b[$j]\n";
}

for ($j=0;$j<$NW;$j++) {
#$l[0]=0;
    $r[$j]=1;
    $g[$j]=1;
    $b[$j]=1;

#    $r[$j]=0.9;
#    $g[$j]=0.9;
#    $b[$j]=0.9;
}
#    $r[0]=1;
#    $g[0]=1;
#    $b[0]=1;
#$r[$nc-1]=1;
#$g[$nc-1]=1;
#$b[$nc-1]=1;

pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

#pgsvp(0.15,0.7,0.1,0.9);
#pgswin($x_min,$x_max,$y_min,$y_max);

#pgenv($x_min,$x_max,$y_min,$y_max,1,0);
pgenv($x_min,$x_max,$y_min,$y_max,1,0);
#pgpap(6.0,1.0);

#pgsch(0.7);           # Set character height
#pgsch(1.2);

if ($dev !~ "TPNG") {
    pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
}
pgimag(\@map,$nx,$ny,1,$nx,1,$ny,$min,$max,\@tr);

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
pgsch(1.3);
#pgptxt(-13,28,0,0.5,$label);
if ($dev !~ "TPNG") {

    if ($rot==0) {
	pgptxt(-15,28,0,0.5,$label);
    } else {
	pgptxt(-15,24,0,0.5,$label);
    }
}
pgsci(1);
pgsch(1);

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
pgslw(3);
if ($dev !~ "TPNG") {
    pgsci(1);
    pgslw(3);
} else {
    if ($reverse==0) {
        pgsci(4);
    } else {
        pgsci(11);
    }
    pgslw(8);
}

pgcons(\@cont,$nx,$ny,1,$nx,1,$ny,\@levels,$nlevels,\@tr);
pgsci(1);
pgslw(1);



for ($i=0;$i<$nlevels2;$i++) {    
    $levels2[$i]=$min_cont2+$steps2*($i**1.5);
#    print "$i $levels[$i]\n";
}

pgsci(12);
pgslw(4);
pgsls(2);
#pgsfs(2);

if ($dev !~ "TPNG") {
    pgsci(1);
    pgslw(3);
} else {
    if ($reverse==0) {
        pgsci(4);
    } else {
        pgsci(3);
    }
    pgslw(8);
}


pgcons(\@cont2,$nx,$ny,1,$nx,1,$ny,\@levels2,$nlevels2,\@tr);
#pgsls(3);
pgsls(1);
pgslw(1);
pgsci(1);

if ($dev !~ "TPNG") {

#pgptxt(-85,160,0,0,$label);
#pgptxt(-110,160,0,0,$label);
pgwedg("RI",0.0,4.5,$min,$max,"");
#pgwedg("RI",0,1.0,1,500,"");
pgbox("SBC",0,0,"SBC",0,0);
}




pgclos();
pgend();



exit;
