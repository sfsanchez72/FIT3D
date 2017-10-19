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
    print "USE: spec2D_plot.pl INPUT_FILE.FITS NY DEVICE [MIN MAX] [WMIN WMAX]\n [REF_LINES CUT]";
   exit;
}

$input=$ARGV[0];
$NY=$ARGV[1];
$dev=$ARGV[2];
$def=0;
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}


@SN=split(/\./,$input);
$SN_input=$SN[0].".SN.fits";




$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;
open(OBJ,"dump_header.pl $input | ");
while($l_obj=<OBJ>) {
    chop($l_obj); 
    ($key,$val)=split(/\|/,$l_obj);
    $key =~ s/ //g;
    $val =~ s/ //g;
    $h{$key}=$val;
}
close(OBJ);


#$h=$pdl->gethdr;
$pdl_SN=rfits($SN_input);
#print "$crval $cdelt $crpix\n";
if ($cdelt==0) {
    $cdelt=1;
}
$k=$NY+1;
$val="VAL1_".$k;
$crval=$pdl->hdr->{$val};
$val="DLT1_".$k;
$cdelt=$pdl->hdr->{$val};
#print "$crval $cdelt\n";

for ($i=0;$i<$nx;$i++) {

#    $crpix=$pdl->hdr->{"CRPIX1"};


    $wave[$i]=$crval+$cdelt*$i;
    $flux[$i]=$pdl->at($i,$NY);
    $SN_flux[$i]=$pdl_SN->at($i,$NY);
    if ($flux[$i]>$y_max) {
	$y_max=$flux[$i];
    }
    if ($flux[$i]<$y_min) {
	$y_min=$flux[$i];
    }
#    print "$crval $cdelt $crpix $wave[$i] $flux[$i]\n";
}
$med_flux=median(@flux);
$y_min=-0.2*$med_flux;
$y_max=3*$med_flux;


$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}

if ($#ARGV==6) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $def=1;
}
$ref_def=0;
if ($#ARGV==8) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $ref_lines=$ARGV[7];
    $cut_lines=$ARGV[8];
    $ref_def=1;
    $def=1;
}

if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}

#print "$y_min $y_max\n";

if ($ref_def==1) {
    $nl=0;
    open(FH,"<$ref_lines");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($data[0]>$cut_lines) {
	    $wave_line[$nl]=$data[1];
	    if ($nl>0) {
		if (abs($wave_line[$nl]-$wave_line[$nl-1])>5) {
		    $nl++;
		    }
	    } else {
		$nl++;
	    }
	}
    }
    close(FH);
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
#pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pgsvp(0.15,0.85,0.35,0.9);
pgswin($wmin,$wmax,$y_min,$y_max);
pglabel("","Flux (units)","");
pgline($nx,\@wave,\@flux);
pgsch(1.0);
pgbox("TSBC",0,0,"SBCNV",0,0);

pgsci(2);
for ($j=0;$j<$nl;$j++) {
    if (($wave_line[$j]>$wmin)&&($wave_line[$j]<$wmax)) {
	if ($j==(2*int($j/2))) {
	    pgptxt($wave_line[$j],0.8*$y_max,90,0,"$wave_line[$j]");
	} else {
	    pgptxt($wave_line[$j],0.7*$y_max,90,0,"$wave_line[$j]");
	}
    }
}
#
# AB mag
#

my $dA=$wave[1]-$wave[0];

$med_w=0.5*($wmax+$wmin);
$med_F=median(@flux);
#$med_F=$med_F*$dA;
$area=1;
$mag=8.9-2.5*log10(($med_F*1e-16*(($med_w)**2)/((1e-23)*(3e18)))/$area);
$med_SN=median(@SN_flux);
$exp=$h{'EXPTIME'};
$obj=$h{'OBJECT'};
$date=$h{'DATE'};
$ord=$NY+60;
pgptxt($wmin+0.05*($wmax-$wmin),$y_max-0.07*($y_max-$y_min),0,0,"'$obj' #$ord");
pgptxt($wmin+0.05*($wmax-$wmin),$y_max-0.14*($y_max-$y_min),0,0,"T\\dexp\\u $exp s");
$cut=apr($mag);
pgptxt($wmin+0.05*($wmax-$wmin),$y_max-0.21*($y_max-$y_min),0,0,"mag\\dAB\\u $cut");
$cut=apr($med_SN);
pgptxt($wmin+0.05*($wmax-$wmin),$y_max-0.28*($y_max-$y_min),0,0,"S/N $cut");

print "$date | $obj | $exp | $med_w | $mag | $med_SN \n";

pgsci(1);
pgsch(1.4);           # Set character height
pgsvp(0.15,0.85,0.15,0.35);
pgswin($wmin,$wmax,-0.5*$med_SN,3*$med_SN);
pglabel("Wavelength (\\A)","S/N","");
pgline($nx,\@wave,\@SN_flux);
pgsch(1.0);
pgbox("NTSBC",0,0,"NSBCV",0,0);


pgclose;
pgend;


exit;
