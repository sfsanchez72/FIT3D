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


if ($#ARGV<2) {
    print "USE: poly.pl INPUT1.FITS NPOLY OUTPUT.FITS [NX1 NX2] [DEV] [LIMIT]\n";
    exit;
}

$infile=$ARGV[0];
$npoly=$ARGV[1];
$outfile=$ARGV[2];
$a_in=rfits($infile);
$h=$a_in->gethdr;
($nx,$ny)=$a_in->dims;
$nx1=0;
$nx2=0;
$dev="/xs";
if ($#ARGV>2) {
    $nx1=$ARGV[3];
    $nx2=$ARGV[4];
    if ($#ARGV>3) {
	$dev=$ARGV[5];
	$limit=$ARGV[6];
    }
} else {
    $nx2=$nx-1;
}
 

$a_out=zeroes($nx,$ny);
for ($j=0;$j<$ny;$j++) {
    my @a_spec;
    my @w_out;
    my @w_in;
    $nt=0;
    $min=1e12;
    $max=-1e12;

    for ($i=0;$i<$nx;$i++) {
	$w_out[$i]=$i;
	$val=$a_in->at($i,$j);
	if (($i>$nx1)&&($i<$nx2)&&($val>$limit)) {
	    $w_in[$nt]=$i;
	    $a_spec[$nt]=$val;
	    
	    if ($min>$a_spec[$nt]) {
		$min=$a_spec[$nt];
	    }
	    if ($max<$a_spec[$nt]) {
		$max=$a_spec[$nt];
	    }
	    
	    $nt++;
	}	 
	
    }

#    print "$nt\n";
    if ($nt>($npoly+10)) {
	($s_x,$coeff) = fitpoly1d(pdl(@w_in),pdl(@a_spec),$npoly);
	@c=list($coeff);
	for ($i=0;$i<$nx;$i++) {
	    $val=$a_in->at($i,$j);
	    if (($i<=$nx1)&&($i=>$nx2)&&($val<=$limit)) {
		$out_spec[$i]=$val;
	    } else {
		$out_spec[$i]=0;
		for ($k=0;$k<$npoly;$k++) {
		    $out_spec[$i]=$out_spec[$i]+$c[$k]*($i**$k);
#	    print "$k $c[$k] $out_spec[$i]\n";
		}

	    }
	    set($a_out,$i,$j,$out_spec[$i]);
	}
    } else {
	for ($i=0;$i<$nx;$i++) {
	    $val=$a_in->at($i,$j);
	    $out_spec[$i]=$val;
	    set($a_out,$i,$j,$out_spec[$i]);
	}
    }

    pgbegin(0,$dev,1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv(0,$nx,$min,$max,0,0);
    pglabel("X-axis","Y-Axis","$j/$ny");
    pgsci(1);
    pgpoint($nt,\@w_in,\@a_spec,5);
    pgsci(3);
    pgline($nx,\@w_out,\@out_spec);
    pgsci(1);
    pgsci(1);
    pgclose;
    if ($dev !~ "null") {
	print "Press Enter\n";
	<stdin>;
    }
}

$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
