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
    print "USE: fiber_trans_flat.pl SKYFLAT.fits FIBERFLAT.fits [NX_min] [NX_max]\n";
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
if ($#ARGV==3) {
    $nx_min=$ARGV[2];
    $nx_max=$ARGV[3];
} else {
    $nx_min=0;
    $nx_max=$nx;
}


my @med;
my @medos;
$k=0;
for ($j=0;$j<$ny;$j++) {
    my @spec;
    for ($i=$nx_min;$i<$nx_max;$i++) {
	$spec[$i]=$in_array[$j][$i];
#	print "$i $spec[$i]\n";
    }
    $med[$j]=mean(@spec);
    if (($med[$j]>0)&&($med[$j]<100000)) {
	$medos[$k]=$med[$j];
	$k++;
    }
#   print "med[$j]=$med[$j]\n";
}
$mean=median(@medos);
#print "med=$mean\n";

for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	if ($mean!=0) {
	    if (($mean>0)&&($med[$j]<100000)) {
		$out_array[$j][$i]=$med[$j]/$mean;
	    } else {
		$out_array[$j][$i]=1/100;
	    }
	} else {
	    $out_array[$j][$i]=0;
	}
    }
}

print "Writting the fiber transmision flat flat $out_file\n";
system("rm $out_file");
write_img($out_file,$nx,$ny,\@out_array);
print "DONE\n";



exit;
