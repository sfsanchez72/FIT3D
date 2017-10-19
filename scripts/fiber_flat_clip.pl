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


if ($#ARGV<1) {
    print "USE: fiber_flat.pl SKYFLAT.fits FIBERFLAT.fits\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$x_smooth=10;
#if ($#ARGV==2) {
#    $x_smooth=$ARGV[2];
#}
#$fiber_flat=$ARGV[2];


print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";

#$smoothed=med2df(pdl(@in_array),$x_smooth,1,{Boundary => Reflect});


my @spec;
for ($i=0;$i<$nx;$i++) {
    my @cut;
    $k=0;
    for ($j=0;$j<$ny;$j++) {
	if ($in_array[$j][$i]>0) {
	    $cut[$k]=$in_array[$j][$i];
#	$cut[$j]=$smoothed->at($i,$j);
	    $k++;
	}
    }
#    $med=median(pdl(@cut));
#    $spec[$i]=median(@cut);

    $spec_tmp=median(@cut);

    my @new_cut;
    $nc=0;
    for ($j=0;$j<$k;$j++) {
	if (abs($cut[$j]-$spec_tmp)<0.15*$spec_tmp) {
	    $new_cut[$nc]=$cut[$j];
	    $nc++;
	}
    }


    $spec[$i]=median(@new_cut);



#    $spec[$i]=median(pdl(@cut));
    for ($j=0;$j<$ny;$j++) {
#	$a=$smoothed->at($i,$j);
#	$out_array[$j][$i]=$a/$med;
	if ($spec[$i]!=0) {
	    $out_array[$j][$i]=$in_array[$j][$i]/$spec[$i];
	} else {
	    $out_array[$j][$i]=1;
	}
    }
}

print "Writting the fiber flat $out_file\n";
system("rm $out_file");
write_img($out_file,$nx,$ny,\@out_array);
print "DONE\n";

exit;
