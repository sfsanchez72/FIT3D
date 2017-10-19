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
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<4) {
    print "USE: cr_lap_clean.pl input.fits nsigma limit clean.fits cr.fits [AXIS=0/1]\n";
    exit;
}

$infile=$ARGV[0];
$nsigma=$ARGV[1];
$limit=$ARGV[2];
$outfile=$ARGV[3];
$mask=$ARGV[4];
$axis=0;
if ($#ARGV==5) {
    $axis=$ARGV[5];
}



$pdl=rfits($infile);
#print "$pdl\n";
$h=$pdl->gethdr;
#$pdl=$pdl*1.01;
#$pdl=$pdl/1.01;
($nx,$ny)=$pdl->dims();

$kernel=pdl([[0,-1,0],[-1,4,-1],[0,-1,0]]);
$smoothed1=med2df($pdl,3,3,{Boundary => Reflect});
#$smoothed2=med2df($pdl,5,5,{Boundary => Reflect});
$s=$smoothed1;
#$s->wfits("smooth.fits");
$l=conv2d($pdl,$kernel,{Boundary => Reflect});
$F=$l/$s;
$F->inplace->setnantobad;
$F->inplace->setbadtoval(0); 


#$F->wfits("F.fits");

$nx1=int($nx/2-200);
$nx2=int($nx/2+200);
$ny1=int($ny/2-200);
$ny2=int($ny/2+200);

if ($nx1<0) {
    $nx1=0;
    $nx2=101;
}

if ($ny1<0) {
    $ny1=0;
    $ny2=101;
}

if ($ny2>$ny-1) {
    $ny2=$ny-1;
    $ny1=$ny-101;
}

if ($nx2>$nx-1) {
    $nx2=$nx-1;
    $nx1=$nx-101;
}

$sec=$F->slice("$nx1:$nx2,$ny1:$ny2");
($mean,$rms,$median,$min,$max)=stats($sec);
#print "$mean,$rms,$median,$min,$max\n";
$pdl_mask=ones($nx,$ny);
$pdl_out=$pdl;
for ($j=1;$j<$ny-1;$j++) {
    for ($i=1;$i<$nx-1;$i++) {	
	$val=$F->at($i,$j);
	$val_pdl=$pdl->at($i,$j);
	if (($val>$nsigma*abs($rms))&&($val>$limit)&&($val_pdl>$limit)) {
	    for ($ii=-1;$ii<2;$ii++) {
		for ($jj=-1;$jj<2;$jj++) {
		    set($pdl_mask,$i+$ii,$j+$jj,0);
		}
	    }
	}
    }
}


if ($axis==0) {
    for ($j=1;$j<$ny-1;$j++) {
	for ($i=1;$i<$nx-1;$i++) {
	    $val_mask=$pdl_mask->at($i,$j);
	    if ($val_mask==0) {
		$last=$pdl_out->at($i-1,$j);
		set($pdl_out,$i,$j,$last);
	    }
	}
    }
} else {
    for ($i=1;$i<$nx-1;$i++) {	
	for ($j=1;$j<$ny-1;$j++) {
	    $val_mask=$pdl_mask->at($i,$j);
	    if ($val_mask==0) {
		$last=$pdl_out->at($i,$j-1);
		set($pdl_out,$i,$j,$last);
	    }
	}
    }

}


#$pdl_out=$F;




$pdl_out->sethdr($h);
$pdl_out->wfits($outfile);
$pdl_mask->sethdr($h);
$pdl_mask->wfits($mask);
$F->wfits("F.fits");


exit;

