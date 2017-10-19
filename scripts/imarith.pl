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


if ($#ARGV<3) {
    print "USE: imarith.pl INPUT1.FITS OPERATOR(+,-,/,*) INPUT2.FITS OUTPUT.FITS\n";
    exit;
}

$infile1=$ARGV[0];
$operator=$ARGV[1];
$infile2=$ARGV[2];
$outfile=$ARGV[3];

$a_in1=rfits($infile1);
$h=$a_in1->gethdr;
if ($infile2 !~ /fit/) {
#$val=$infile2*1.0;
#print "$infile2 $val\n";
#if (($infile2 == $val)&&($val!=0)) {  
    $a_in2=$infile2;
} else {
    $a_in2=rfits($infile2);
}

if ($operator eq "+") {
    $a_in1=$a_in1+$a_in2;
#    print "+\n";
}

if ($operator eq "-") {
    $a_in1=$a_in1-$a_in2;
#    print "-\n";
}

if ($operator eq "/") {
    $a_in1=$a_in1/$a_in2;
#    print "/\n";
}

#print "$operator\n";
if ($operator eq "*") {
    $a_in1=$a_in1*$a_in2;
}
$a_in1->sethdr($h);
$a_in1->wfits($outfile);

exit;
