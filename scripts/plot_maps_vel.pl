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

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<7) {
    print "USE: plot_maps.pl map.fits min max bright contrast Label factor dev [MASK_MAP]\n";
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
$id_mask=0;
if ($#ARGV==8) {
    $mask_file=$ARGV[8];
    $id_mask=1;
}


#$table="idl5"; $reverse=1; 
$table="smooth2"; $reverse=0; 
#$table="real"; $reverse=0; 
#$bright=0.7; 
#$contra=0.5; 
$color_cont=0;

#$pdl_mask=rfits("M74_mos.mask.fits");


$pdl_map=rfits($mapfile);
#$crval1=$pdl_map->hdr->{CRVAL1};
#$cdelt1=$pdl_map->hdr->{CDELT1};
#$crval2=$pdl_map->hdr->{CRVAL2};
#$cdelt2=$pdl_map->hdr->{CDELT2};
#print "$crval1 $cdelt1 $crval2 $cdelt2\n";
($nx,$ny)=$pdl_map->dims;
$h=$pdl_map->gethdr;
$crval1=-138.839996337891;
$cdelt1=4;
$crval2=-210.726806640625;
$cdelt2=4;

if ($id_mask==1) {
    $pdl_mask=rfits($mask_file);
    $pdl_map=$pdl_map*$pdl_mask*$factor;
} else {
    $pdl_map=$pdl_map*$factor;
}


$NNx=$nx-1;
$NNy=$ny-1;
$rpdl_map=$pdl_map->slice("$NNx:0,$NNy:0");


@map=list($rpdl_map);
$crval1=$crval1-50-10;
$crval2=$crval2+17-10;
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
for ($j=0;$j<128;$j++) {
    $l[$j] = $j*(1/256);
    $r[$j] = 1;
    $g[$j] = $j*(2/256);
    $b[$j] = $j*(2/256);
}

for ($j=128;$j<256;$j++) {
    $l[$j] = $j*(1/256);
    $r[$j] = (1-$b[$j-128]);
    $g[$j] = (1-$g[$j-128]);
    $b[$j] = 1;
}


pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);


pgenv($x_min,$x_max,$y_min,$y_max,1,0);
#pgpap(6.0,1.0);

#pgsch(0.7);           # Set character height
#pgsch(1.2);

pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","$label");
pgimag(\@map,$nx,$ny,1,$nx,1,$ny,$min,$max,\@tr);

pgsci(5);
pgsch(1.2);
#pgsch(1.1);
@tmp;
$nt=0;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$map[$i+$j*$nx];
#	if ($map[$i+$j*$nx]<=0) {
	if ($val==0) {
	    $X=$crval1+$cdelt1*$i+2;
	    $Y=$crval2+$cdelt2*$j+2;
	#    print "$i,$j $X,$Y MAP=$map[$i+$j*$nx]\n";
	    pgpoint(1,[$X],[$Y],16);
	} else {
#	    if ($val>1) {
		$tmp[$nt]=$val;
#		print "$tmp[$nt] $nt\n";
		$nt++;
#	    }
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
#pgpoint(1,[3.65],[-3.13],3);
pgpoint(1,0,0,3);
pgsch(1.2);
pgslw(1);

pgsci(1);
#pgptxt(-85,160,0,0,$label);
#pgptxt(-110,160,0,0,$label);
#pgptxt(-100,160,0,0,$label);
pgsci(1);
pgwedg("RI",0.0,4.5,$min,$max,"");
#pgwedg("RI",0,1.0,1,500,"");
pgbox("SBC",0,0,"SBC",0,0);





pgclos();
pgend();



exit;
