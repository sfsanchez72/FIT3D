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
    print "USE: poly_spec_clip.pl input.spec NPOLY OUTPUT.spec [NX1 NX2] [DEV]\n";
    exit;
}

$infile=$ARGV[0];
$npoly=$ARGV[1];
$outfile=$ARGV[2];

$nx=0;
open(FH,"<$infile");
while($line=<FH>) {
    @data=split(" ",$line);
    if ($line !~ "#") {
	$index_unc[$nx]=$data[0];
	$wave_unc[$nx]=$data[1];
	$flux_unc_ini[$nx]=$data[2];
	$nx++;
    }
}
close(FH);

$med=median(@flux_unc_ini);
$sig=sigma(@flux_unc_ini);
for ($i=0;$i<$nx;$i++) {
    if (abs($flux_unc_ini[$i]-$med)>0.5*$sig) {
	$flux_unc_ini[$i]=$med;
    }
}


$nx1=0;
$nx2=0;
$dev="/xs";
if ($#ARGV>2) {
    $nx1=$ARGV[3];
    $nx2=$ARGV[4];
    if ($#ARGV>3) {
	$dev=$ARGV[5];
    }
} else {
    $nx2=$nx-1;
}
 
$a_out=zeroes($nx);
my @a_spec;
my $w_out;
$nt=0;
$min=1e12;
$max=-1e12;

for ($i=0;$i<$nx;$i++) {
    $w_out[$i]=$i;
    if (($i>$nx1)&&($i<$nx2)) {
	$w_in[$nt]=$i;
	$a_spec[$nt]=$flux_unc_ini[$i];
	
	if ($min>$a_spec[$nt]) {
	    $min=$a_spec[$nt];
	}
	if ($max<$a_spec[$nt]) {
	    $max=$a_spec[$nt];
	}
	
	$nt++;
    }	 
    
}
($s_x,$coeff) = fitpoly1d(pdl(@w_in),pdl(@a_spec),$npoly);
@c=list($coeff);
for ($i=0;$i<$nx;$i++) {
    $out_spec[$i]=0;
    for ($k=0;$k<$npoly;$k++) {
	    $out_spec[$i]=$out_spec[$i]+$c[$k]*($i**$k);
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


open(FH,">$outfile");
for ($i=0;$i<$nx;$i++) {
    print FH "$index_unc[$i] $wave_unc[$i] $out_spec[$i]\n";
}
close(FH);

exit;
