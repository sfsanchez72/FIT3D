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
    print "USE: plot_templates.pl INPUT_FILE.FITS NY DEVICE [MIN MAX] [WMIN WMAX]\n [REF_LINES CUT]";
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


$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL1"};
$cdelt=$pdl->hdr->{"CDELT1"};
$crpix=$pdl->hdr->{"CRPIX1"};
if ($cdelt==0) {
    $cdelt=1;
}
for ($i=0;$i<$nx;$i++) {
    $wave[$i]=$crval+$cdelt*($i+1-$crpix);
    $flux[$i]=$pdl->at($i,$NY);
    if ($flux[$i]>$y_max) {
	$y_max=$flux[$i];
    }
    if ($flux[$i]<$y_min) {
	$y_min=$flux[$i];
    }
#    print "$wave[$i] $flux[$i]\n";
}




for ($iii=0;$iii<$ny;$iii++) {
    $header="NAME".$iii;
    $name[$iii]=$pdl->hdr->{$header};
    $name_min=$name[$iii];
    $name_min =~ s/spec_ssp_//;
    $name_min =~ s/.spec//;
    ($AGE,$MET)=split("_",$name_min);
    if ($AGE =~ "Myr") {
	$age=$AGE;
	$age =~ s/Myr//;
	$age=$age/1000;
    } else {
	$age=$AGE;
	$age =~ s/Gyr//;
    }
    $met=$MET;
    $met =~ s/z/0\./;

    $age_mod[$iii]=$age;
    $met_mod[$iii]=$met;
    $header="NORM".$iii;
 $val_ml=$pdl->hdr->{$header};
    if ($val_ml!=0) {
	$ml[$iii]=1/$val_ml;
    } else {
	$ml[$iii]=1;
    }
}



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
pgsch(1.2);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
pglabel("Wavelength","Flux","");
$color=1;
for ($NY=0;$NY<$ny;$NY++) {
    for ($i=0;$i<$nx;$i++) {
	$flux[$i]=$pdl->at($i,$NY);
    }
    pgsci($color);
    pgline($nx,\@wave,\@flux);
    $color++;
    if ($color>15) {
	$color=1;
    }
}
pgsch(0.9);
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
pgsci(0);
$dw=$wmax-$wmin;
pgrect($wmin+0.7*$dw,$wmin+0.99*$dw,0.7*$y_max,0.97*$y_max);
pgsci(1);
pgsfs(2);
pgrect($wmin+0.7*$dw,$wmin+0.99*$dw,0.7*$y_max,0.97*$y_max);

pgslw(3);
pgsch(0.9);
$color=1;
for ($NY=0;$NY<$ny;$NY++) {
    $F=(0.72+0.25*$NY/$ny)*$y_max;#*($y_max-$y_min);
    $X1=$wmin+0.71*($wmax-$wmin);
    $X2=$wmin+0.91*($wmax-$wmin);
    $X3=$wmin+0.93*($wmax-$wmin);

    pgsci($color);
#    pgline(2,[$X1,$X2],[$F,$F]);
#    print "$X1 $X2 $F\n";
    $call=$age_mod[$NY]." Gyr ".$met_mod[$NY]." met";
    pgptxt($X1,$F,0,0,$call);
    $color++;
    if ($color>15) {
	$color=1;
    }
}



#pgbox();



pgclose;
pgend;

exit;
