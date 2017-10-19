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
    print "USE: fiber_flat.pl SKYFLAT1.fits SKYFLAT2.fits FIBERFLAT1.fits FIBERFLAT2.fits [NX_min] [NX_max]\n";
    exit;
}

$infile1=$ARGV[0];
$infile2=$ARGV[1];
$out_file1=$ARGV[2];
$out_file2=$ARGV[3];

#$fiber_flat=$ARGV[2];


print "Reading file $infile1\n ";
$nax=read_naxes($infile1);   
@naxis=@$nax;
$nx1=$naxis[0];
$ny1=$naxis[1];
@in_array1=read_img($infile1);
print "DONE\n";
if ($#ARGV==3) {
    $nx_min=$ARGV[4];
    $nx_max=$ARGV[5];
} else {
    $nx_min=0;
    $nx_max=$nx1;
}

print "Reading file $infile2\n";
$nax=read_naxes($infile2);   
@naxis=@$nax;
$nx2=$naxis[0];
$ny2=$naxis[1];
@in_array2=read_img($infile2);
if ($nx1!=$nx2) {
    print "FILES $infile1 and $infile2 do not cover the same range: $nx1!=$nx2\n";
    exit;
}
$nx=$nx1;
print "DONE\n";

my @med;
for ($j=0;$j<$ny1;$j++) {
    my @spec;
    for ($i=$nx_min;$i<$nx_max;$i++) {
	$spec[$i]=$in_array1[$j][$i];
    }
    $med[$j]=median(@spec);
}
for ($j=$ny1;$j<($ny1+$ny2);$j++) {
    my @spec;
    for ($i=$nx_min;$i<$nx_max;$i++) {
	$spec[$i]=$in_array2[$j-$ny1][$i];
    }
    $med[$j]=median(@spec);
}
$mean=mean(@med);

for ($j=0;$j<$ny1;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$out_array1[$j][$i]=$med[$j]/$mean;
    }
}
for ($j=$ny1;$j<($ny1+$ny2);$j++) {
    for ($i=0;$i<$nx;$i++) {
	$out_array2[$j-$ny1][$i]=$med[$j]/$mean;
    }
}


print "Writting the fiber transmision flats $out_file1 & $out_file2\n";
system("rm $out_file1");
write_img($out_file1,$nx,$ny1,\@out_array1);
system("rm $out_file2");
write_img($out_file2,$nx,$ny2,\@out_array2);
print "DONE\n";

exit;
