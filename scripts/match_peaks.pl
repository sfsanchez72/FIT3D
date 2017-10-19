#!/usr/bin/perl
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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<9) {
    print "USE: match_peaks.pl REFERENCE.fits INPUT.fits Spectral_axis[0/1] plot i_center i_width j_center j_width y_shift_init y_shift_limit\n";
    exit;
}

$reffile=$ARGV[0];
$infile=$ARGV[1];
$spec_axis=$ARGV[2];
$plot=$ARGV[3];
$nspec=$ARGV[4];
$width=$ARGV[5];
$j_center=$ARGV[6];
$j_width=$ARGV[7];
$y_shift_init=$ARGV[8];
$y_shift_limit=abs($ARGV[9]);
$ref_pdl=rfits($reffile);
($ref_nx,$ref_ny)=$ref_pdl->dims();
$in_pdl=rfits($infile);
($nx,$ny)=$in_pdl->dims();

if (($nx!=$ref_nx)&&($ny!=$ref_ny)) {
    print "Reference file and input file do not the same dimensions\n";
    exit;
}

$j_min=$j_center-$j_width;
$j_max=$j_center+$j_width;
$i_min=$nspec-$width;
$i_max=$nspec+$width;

if ($spec_axis==0) {
    if ($j_min<0) {
	$j_min=0;
    }
    if ($j_max>=$ny) {
	$j_max=$ny-1;
    }
    if ($i_min<0) {
	$i_min=0;
    }
    if ($i_max>=$nx) {
	$i_max=$nx-1;
    }
    $tmp_ref=$ref_pdl->slice("$i_min:$i_max,$j_min:$j_max");
    $jj_min=$j_min+$y_shift_init;
    $jj_max=$j_max+$y_shift_init;
    $tmp_in=$in_pdl->slice("$i_min:$i_max,$jj_min:$jj_max");
    $sec_ref=sumover($tmp_ref);
    $sec_in=sumover($tmp_in);    
#    $sec_ref=$tmp2_ref->xchg(0,1);
#    $sec_in=$tmp2_in->xchg(0,1);
} else {
    if ($j_min<0) {
	$j_min=0;
    }
    if ($j_max>=$nx) {
	$j_nax=$nx-1;
    }
    if ($i_min<0) {
	$i_min=0;
    }
    if ($i_max>=$ny) {
	$i_max=$ny-1;
    }

    $tmp_ref=$ref_pdl->slice("$j_min:$j_max,$i_min:$i_max");
    $jj_min=$j_min+$y_shift_init;
    $jj_max=$j_max+$y_shift_init;
    $tmp_in=$in_pdl->slice("$jj_min:$jj_max,$i_min:$i_max");
    $tmp2_ref=sumover($tmp_ref->xchg(0,1));
    $tmp2_in=sumover($tmp_in->xchg(0,1));
    $sec_ref=$tmp_ref;#->xchg(0,1);
    $sec_in=$tmp_in;#->xchg(0,1);

}
@box=($j_min,$j_max,$min,$max);
@x=($j_min .. $j_max);
if ($plot==1) {
    @a_ref=list($sec_ref);
    @a_in=list($sec_in);
    ($tmin,$max)=minmax(@a_ref);
    ($min,$tmax)=minmax(@a_in);
    $rat=$max/$tmax;
    @a_in=list($sec_in*$rat);
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgenv($j_min,$j_max,$min,$max,0,0);
    pglabel("Y-axis","Counts","");
    pgline($j_max-$j_min,\@x,\@a_ref);    
    pgsci(3);
    pgline($j_max-$j_min,\@x,\@a_in);    
    pgsci(1);
    pgclose;
    pgend;
}

exit;
