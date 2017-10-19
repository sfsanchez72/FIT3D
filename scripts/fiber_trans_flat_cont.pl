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
    print "USE: fiber_trans_flat.pl SKYFLAT.fits FIBERFLAT.fits NX_min_E NX_max_E NX_min_con NX_max_con\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$nx_min=$ARGV[2];
$nx_max=$ARGV[3];
$nxc_min=$ARGV[4];
$nxc_max=$ARGV[5];



#$fiber_flat=$ARGV[2];


print "Reading file $infile ";
$nax=read_naxes($infile);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
@in_array=read_img($infile);
print "DONE\n";




my @med;
my @medos;
$k=0;
for ($j=0;$j<$ny;$j++) {
    my @spec;
    my @cont;
    for ($i=$nx_min;$i<$nx_max;$i++) {
	$spec[$i]=$in_array[$j][$i];
    }
    for ($i=$nxc_min;$i<$nxc_max;$i++) {
	$cont[$i]=$in_array[$j][$i];
    }
    $med=mean(@spec);
    $medc=mean(@cont);
    $MM[$j]=$med-$medc;
    if (($MM[$j]>0)&&($MM[$j]<1e12)) {
	$medos[$k]=$MM[$j];
	$k++;
    }
}
$mean=median(@medos);
#print "med=$mean\n";

for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	if ($mean!=0) {
	    if (($mean>0)&&($MM[$j]<1e12)) {
		$out_array[$j][$i]=$MM[$j]/$mean;
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
