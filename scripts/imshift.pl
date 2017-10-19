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
use PDL::Transform;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<3) {
    print "USE: imshift.pl INPUT.FITS OUTPUT.FITS DX DY\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$dx=$ARGV[2];
$dy=$ARGV[3];
$shift=pdl($dx,$dy);

$a_in=rfits($infile);
$h=$a_in->gethdr;
#$tf = t_fits($a_in);
#$tr = t_offset($dx,$dy);
$tr = t_offset($shift);
$a_out=$a_in->map($tr,{pix=>1});
$a_out->sethdr($h);
$a_out->wfits($outfile);

exit;
