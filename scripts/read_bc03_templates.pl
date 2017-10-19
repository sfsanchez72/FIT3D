#!/usr/bin/perl
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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<5) {
    print "USE read_bc03_templates.pl template CRVAL CDELT NPIX output_spec.txt MEDIAN_BOX\n";
    exit;
}

$infile=$ARGV[0];
$crval=$ARGV[1];
$cdelt=$ARGV[2];
$npix=$ARGV[3];
$outfile=$ARGV[4];
$med_box=$ARGV[5];

$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    if ($line !~ "#") {
	@data=split(" ",$line);
	$wave_in[$n]=$data[0];
	$flux_in[$n]=$data[1];
	$n++;
    }
}
close(FH);

for ($j=0;$j<$npix;$j++) {
    $wave_out[$j]=$crval+$cdelt*$j;    
}


my $flux_out = interpol(pdl(@wave_out), pdl(@wave_in), pdl(@flux_in));

if ($med_box>1) {
    $a=ones($med_box);
    $flux_out=$flux_out->conv1d($a);
}


open(OUT,">$outfile");
for ($j=0;$j<$npix;$j++) {
    $k=$j+1;
    $val=$flux_out->at($j);
    print OUT "$k $wave_out[$j] $val\n";
}
close(OUT);

exit;
