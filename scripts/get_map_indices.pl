#!/usr/bin/perl
#
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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<0) {
    print "USE: get_map_indices.pl INDICES.CUBE.fits\n";
    exit;
}

$input_cube=$ARGV[0];
$pre=$input_cube;
$pre =~ s/.cube.fits/.fits/;

$b=rfits("$input_cube");
$h=$b->gethdr;
($nx,$ny,$nz)=$b->dims;

for ($k=0;$k<$nz;$k++) {
    my $out=$b->slice(":,:,($k)");
    my $NAME="INDEX_".$k;
    my $head=$b->hdr->{$NAME};
    my $outfile=$head.".".$pre;
    $out->wfits($outfile);    
}

exit;





