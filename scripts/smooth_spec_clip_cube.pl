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

if ($#ARGV<5) {
    print "USE: smooth_spec_clip_cube.pl intput_spec.cube.fits output_spec.cube.fits WIDTH NSIGMA NXMIN NXMAX\n";
    exit;
}

$input=$ARGV[0];
$outfile=$ARGV[1];
$box=$ARGV[2];
$nsigma=$ARGV[3];
$nxmin=$ARGV[4];
$nxmax=$ARGV[5];


$pdl_in=rfits($input);
($nx,$ny,$nz)=$pdl_in->dims;

for ($i=0;$i<$nx;$i++) {
    my $pdl_now=$pdl_in->slice("($i),,");
    @N=$pdl_now->dims;
    my $pdl_cut=$pdl_now->slice(",$nxmin:$nxmax");
    my @stats=stats($pdl_cut);
#    for ($j=0;$j<$ny;$j++) {
#	my $pdl_now=$pdl_in->slice("($i),($j),");
#	my @flux=list($pdl_now);
#	my @Fcut=list($pdl_now->slice("$nxmin:$nxmax"));
    my $med=$stats[2];
    my $sig=$stats[1];
    my $min=$stats[3];
    my $max=$med+$nsigma*$sig;
    $pdl_now=$pdl_now->clip($min,$max);
    my $pdl_smooth=med2df($pdl_now,1,$box);
 #   my $npoly=11;
    
 #   for ($j=0;$j<$ny;$j++) {
#	my $pdl_data=$pdl_smooth->slice("($j),");
#	my $s_x = fitpoly1d($pdl_data,$npoly);
#	$pdl_data .= $s_x;
#    }
    my $pdl_now=$pdl_in->slice("($i),,");   
    $pdl_now .= $pdl_smooth;
#    $pdl_now .= $s_x;
#    print "$pdl_now\n";
#    print "$i/$nx\n";
}

$pdl_in->wfits($outfile);

exit;


