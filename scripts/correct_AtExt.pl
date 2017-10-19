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
    print "USE: correct_AtExt.pl RSS.fits OUT.fits Av Airmass CRVAL CDELT\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$av=$ARGV[2];
$x=$ARGV[3];
$crval=$ARGV[4];
$cdelt=$ARGV[5];

#$fiber_flat=$ARGV[2];


print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";

my @spec;
for ($i=0;$i<$nx;$i++) {
    my @cut;
    $w=$crval+$cdelt*$i;
    $ext_w=0.00533366480746607*($w/10000)**(-4)+0.0807207505773469*($w/10000)**(-0.8)+0.0222704006467989;
    $ext_w=($ext_w/0.23)*$av*$x;
    $factor=10**(0.4*$ext_w);
    for ($j=0;$j<$ny;$j++) {
	$out_array[$j][$i]=$in_array[$j][$i]*$factor;
    }
#    print "$w $factor $ext_w\n";
}

print "Writting the corrected spectra $out_file\n";
system("rm $out_file");
write_rss($out_file,$nx,$ny,$crval,$cdelt,\@out_array);
print "DONE\n";

exit;
