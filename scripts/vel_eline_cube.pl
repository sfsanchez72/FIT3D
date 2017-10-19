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

$vel_light=299792.458;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");

if ($#ARGV<5) {
    print "USE: vel_elines_cube.pl CUBE.fits nsearch LIMIT(% max) PREFIX_OUT WAVE_REF DEV [WMIN WMAX]\n";
    print "OUTPUT: PREFIX_OUT.vel_map.fits\n";
    print "OUTPUT: PREFIX_OUT.mask_map.fits\n";
    exit;
}

$input=$ARGV[0];
$nsearch=$ARGV[1];
$imin=$ARGV[2];
$outfile=$ARGV[3];
$wave_ref=$ARGV[4];
$dev=$ARGV[5];

$map_outfile=$outfile.".vel_map.fits";
$mask_outfile=$outfile.".mask_map.fits";


if ($#ARGV==7) {
    $wmin=$ARGV[6];
    $wmax=$ARGV[7];
} else {
    $wmin=0;
    $wmax=1e12;
}
$dmin=0;


$pdl_in=rfits("$input");
$h=$pdl_in->gethdr;
$crval3=$pdl_in->hdr->{CRVAL3};
$cdelt3=$pdl_in->hdr->{CDELT3};
$crpix3=$pdl_in->hdr->{CRPIX3};
($nx,$ny,$nz)=$pdl_in->dims;
$map_out=zeroes($nx,$ny);
$mask_out=zeroes($nx,$ny);

for ($ii=0;$ii<$nx;$ii++) {
    for ($jj=0;$jj<$ny;$jj++) {

# Define data
	$n=0;
	my @wave;
	my @flux;
	$y_min=1e12;
	$y_max=-1e12;
	for ($k=0;$k<$nz;$k++) {
	    $wave[$n]=$crval3+$cdelt3*$k;
	    $flux[$n]=$pdl_in->at($ii,$jj,$k);
	    if (($wave[$n]>$wmin)&&($wave[$n]<$wmax)) {		
		if ($flux[$n]>$y_max) {
		    $y_max=$flux[$n];
		}
		if ($flux[$n]<$y_min) {
		    $y_min=$flux[$n];
		}
		$n++;
	    }
	}


# Look for peaks
	#$wmin=$wave[0];
	#$wmax=$wave[$n-1];
	$crval=$wave[0];
	$cdelt=$wave[1]-$wave[0];
	$npeaks=0;
	for ($j=$nsearch;$j<($n-$nsearch);$j++) {
	    $peak=1;
	    for ($i=0;$i<$nsearch;$i++) {
		if ($flux[$j-$i]<$flux[$j-$i-1]) {
		    $peak=0;
		}
		if ($flux[$j+$i]<$flux[$j+$i+1]) {
		    $peak=0;
		}
	    }
	    if ($peak==1) {
		$rrr=$imin*$y_max;
		if ($flux[$j]<($imin*$y_max)) {
		    $peak=0;
		}
	    }
	    
	    if ($peak==1) {
		if ($npeaks>0) {
		    $delta=$j-$peak_y_pixel[$npeaks-1];
		    if ($delta<$dmin) {
			$peak=0;
		    }
		}
	    }
	    
	    if ($peak==1) {
		$peak_y_pixel[$npeaks]=$j;
		
		
		$a=$j-1;
		$b=$j;
		$c=$j+1;
		$fa=-$flux[$a];
		$fb=-$flux[$b];
		$fc=-$flux[$c];
		$den=($fc-2*$fb+$fa);
		if ($den!=0) {
		    $peak_y_max[$npeaks]=$c-($b-$a)*(($fc-$fb)/$den+0.5);
		} else {
		    $peak_y_max[$npeaks]=0;
		}
		$npeaks++;
	    }
	    
	}

#
# We plot the section
#
	if (($y_max>$y_min)&&($y_min!=0)) {
	pgbegin(0,"$dev",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($wmin,$wmax,$y_min,$y_max,0,0);
	pglabel("Y-axis","Counts","");
	pgline($n,\@wave,\@flux);    
	pgsch(1.0);           # Set character height
	for ($k=0;$k<$npeaks;$k++) {
	    $wave_peak[$k]=$crval+$cdelt*$peak_y_max[$k];
	    pgsci(2);
	    $x=$crval+$cdelt*$peak_y_pixel[$k];
	    $y=0.8*$y_max;
	    pgpoint(1,[$x],[$y],5);
	    pgsci(4);
	    $x=$crval+$cdelt*$peak_y_max[$k];
	    $y=0.6*$y_max;
	    pgpoint(1,[$x],[$y],2);
	}
	pgsci(3);
	pgline(2,[$wmin,$wmax],[$imin*$y_max,$imin*$y_max]);
	pgsci(1);
	pgsch(1.2);           # Set character height
#    }
	pgsci(1);
	pgclose;
	pgend;
#	print "Press Enter"; <stdin>;
	if ($npeaks==1) {
	    $vel=($wave_peak[0]/$wave_ref-1)*$vel_light;
	    set($map_out,$ii,$jj,$vel);
	    set($mask_out,$ii,$jj,1);
	}
	if ($npeaks==2) {
	    $vel=($wave_peak[0]/$wave_ref-1)*$vel_light;
	    set($map_out,$ii,$jj,$vel);
	    set($mask_out,$ii,$jj,2);
	}
	if ($npeaks==3) {
	    $vel=($wave_peak[1]/$wave_ref-1)*$vel_light;
	    set($map_out,$ii,$jj,$vel);
	    set($mask_out,$ii,$jj,2);
	}
	}
    }
    print "$ii/$nx\n";
    
}

$h = {NAXIS=>2, NAXIS1=>$nx, NAXIS=>$ny, COMMENT=>"vel_eline_cube result"};
$$h{WMIN}=$wmin;
$$h{WMAX}=$wmax;
$$h{WREF}=$wave_ref;
$$h{FILENAME}=$map_outfile;
$map_out->sethdr($h);
$map_out->wfits($map_outfile);
$$h{FILENAME}=$mask_outfile;
$mask_out->sethdr($h);
$mask_out->wfits($mask_outfile);

exit;



