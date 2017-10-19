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
    print "USE: table_stats.pl ASCII_TABLE.txt N.COLUMN(0...N) DEV NBIN [GREP]\n";
    exit;
}

$infile=$ARGV[0];
$nc=$ARGV[1];
$dev=$ARGV[2];
$nbin=$ARGV[3];
$no_grep=0;
if ($#ARGV==4) {
    $no_grep=1;
    $grep=$ARGV[4];
}
$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	if ($no_grep==0) {
	    @data=split(",",$line);
	    $val[$n]=$data[$nc];
	    $n++;
	} else {
	    if ($line =~ $grep) {
		@data=split(",",$line);
		$val[$n]=$data[$nc];
		$n++;
	    }
	}
    }
}
close(FH);
$slice=pdl(@val);
($mean,$rms,$median,$min,$max) = stats($slice);
$sum=sum($slice);
print "$mean $rms $median $min $max $n $sum\n";

$N=$n;
@a=list($slice);
$NN=$n/int($nbin/6+1);


pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
$x_min=$median-3*$rms;
$x_max=$median+3*$rms;
if ($#ARGV==6) {
    $no_grep=1;
    $grep=$ARGV[4];
    $x_min=$ARGV[5];
    $x_max=$ARGV[6];
}
$y_min=0;
$y_max=$NN;
if ($#ARGV==8) {
    $no_grep=1;
    $grep=$ARGV[4];
    $x_min=$ARGV[5];
    $x_max=$ARGV[6];
    $y_min=$ARGV[7];
    $y_max=$ARGV[8];
}
pgenv($x_min,$x_max,$y_min,$y_max,0,0);
pglabel("$infile Value","N values","$median+-$rms");
pgsfs(3);
pghist($N,\@a,$x_min,$x_max,$nbin,3);
pgsfs(2);
pghist($N,\@a,$x_min,$x_max,$nbin,3);
pgclose;
pgend;
exit;





