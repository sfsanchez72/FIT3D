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
    print "USE: sigma_spec.pl RSS.fits SIGMA_SPEC.fits [NY]\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
#$fiber_flat=$ARGV[2];


print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";
$new_ny=$ny;
if ($#ARGV==2) {
    $new_ny=$ARGV[2];
}


my @spec;
for ($i=0;$i<$nx;$i++) {
    my @cut;
    for ($j=0;$j<$ny;$j++) {
	$cut[$j]=$in_array[$j][$i];
    }
    $spec[$i]=sigma(@cut);
    for ($j=0;$j<$new_ny;$j++) {
	$out_array[$j][$i]=abs($spec[$i]);
    }
}

print "Writting the median spectrum $out_file\n";
system("rm $out_file");
write_img($out_file,$nx,$new_ny,\@out_array);
print "DONE\n";

exit;
