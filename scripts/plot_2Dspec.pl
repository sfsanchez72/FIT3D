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
    print "USE: plot_maps.pl map.fits min max bright contrast Label factor dev [TABLE REVERSE] [x_min x_max]\n";
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

$table="idl5"; $reverse=1; 

if ($#ARGV==9) {
     $table=$ARGV[8];
    $reverse=$ARGV[9];
}

#$table="smooth2"; $reverse=0; 
#$table="real"; $reverse=0; 
#$bright=0.7; 
#$contra=0.5; 
$color_cont=0;

#$pdl_mask=rfits("M74_mos.mask.fits");


$pdl_map=rfits($mapfile);
$crval1=$pdl_map->hdr->{CRVAL1};
$cdelt1=$pdl_map->hdr->{CDELT1};
#$crval2=$pdl_map->hdr->{CRVAL2};
#$cdelt2=$pdl_map->hdr->{CDELT2};
#print "$crval1 $cdelt1 $crval2 $cdelt2\n";
($nx,$ny)=$pdl_map->dims;
$h=$pdl_map->gethdr;
#$crval1=-138.839996337891;
#$cdelt1=4;
$crval2=1;
$cdelt2=1;

$pdl_map=$pdl_map*$factor;





$NNx=$nx-1;
$NNy=$ny-1;
#$rpdl_map=$pdl_map->slice("$NNx:0,$NNy:0");


@map=list($pdl_map);
$x_max=$crval1+$cdelt1*$nx;
$x_min=$crval1;
$y_min=1;
$y_max=$ny;


if ($#ARGV==11) {
    $table=$ARGV[8];
    $reverse=$ARGV[9];
    $x_min=$ARGV[10];
    $x_max=$ARGV[11];
}

for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	if ($map[$i+$j*$nx]>1e12) {
	    $map[$i+$j*$nx]=0;
	}
	if ($map[$i+$j*$nx] eq "nan") {
	    $map[$i+$j*$nx]=0;
#
	}
	if ($map[$i+$j*$nx]!=(1*$map[$i+$j*$nx])) {
	    $map[$i+$j*$nx]=0;
	}
    }
}


@tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);




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
}

pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

#pgsvp(0.15,0.7,0.1,0.9);
#pgswin($x_min,$x_max,$y_min,$y_max);

#pgenv($x_min,$x_max,$y_min,$y_max,1,0);
pgenv($x_min,$x_max,$y_min,$y_max,0,0);
#pgpap(6.0,1.0);

#pgsch(0.7);           # Set character height
#pgsch(1.2);

pglabel("Wavelength (\\A)","Spec Id.","$label");
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
