#!/usr/bin/perl
use PDL;
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
use PDL::Func;
use PDL::Math;


if ($#ARGV<3) {
    print "USE: spec_arith.pl INPUT1.spec OPERATOR(+,-,/,*) INPUT2.spec OUTPUT.spec\n";

    exit;
}

$infile1=$ARGV[0];
$operator=$ARGV[1];
$infile2=$ARGV[2];
$output=$ARGV[3];

$n1=0;
open(FH,"<$infile1");
while($line=<FH>) {
   chop($line);
   if ($line !~ "#") {
       @data=split(" ",$line);
       $wave1[$n1]=$data[1];
       $spec1[$n1]=$data[2];
       $n1++;
   }
}
close(FH);

$n2=0;
open(FH,"<$infile2");
while($line=<FH>) {
   chop($line);
   if ($line !~ "#") {
       @data=split(" ",$line);
       $wave2[$n2]=$data[1];
       $spec2[$n2]=$data[2];
       $n2++;
   }
}
close(FH);

if (($wave2[0]!=$wave1[0])||($n1!=$n2)) {
    print "Wavelength range does not match ($wave2[0]!=$wave1[0])\n"; exit;
}

$a_in1=pdl(@spec1);
$a_in2=pdl(@spec2);

if ($operator eq "+") {
    $a_in1=$a_in1+$a_in2;
}

if ($operator eq "-") {
    $a_in1=$a_in1-$a_in2;
}

if ($operator eq "/") {
    $a_in1=$a_in1/$a_in2;
}

if ($operator eq "*") {
    $a_in1=$a_in1*$a_in2;
}

open(OUT,">$output");
#print OUT "# SPEC_ARITH $infile1 $infile2\n";
for ($i=0;$i<$n1;$i++) {
    $val=$a_in1->at($i);
    print OUT "$i $wave1[$i] $val\n";
}
close(OUT);

exit;

