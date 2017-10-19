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

if ($#ARGV<2) {
    print "USE: rebin_cube.pl INPUT_CUBE.fits OUTOUT_CUBE.fits NBIN>1\n";
    exit;
}

$input_cube=$ARGV[0];
$out_file=$ARGV[1];
$nbin=$ARGV[2];

$nbin=int($nbin);
if ($nbin==0) {
    exit;
}

$b=rfits("$input_cube");
$h=$b->gethdr;
$crval3=$b->hdr->{CRVAL3};
$cdelt3=$b->hdr->{CDELT3};
($nx,$ny,$nz)=$b->dims;
if ($nbin>$nz) {
    print "NBIN=$nbin>$nz\n";
    exit;
}
$new_nz=int($nz/$nbin);
$pdl=zeroes($nx,$ny,$new_nz);
$cdelt3=$cdelt3*$nbin;

for ($i=0;$i<$new_nz-1;$i++) {
    $start_i=$i*$nbin;
    $end_i=($i+1)*$nbin;
    my $a=$b->slice(",,$start_i:$end_i");
    my $c=average($a->xchg(0,2));
    my $d=$c->xchg(0,1);
    my $t=$pdl->slice(":,:,$i");
    $t .= $d;
}
$h->{CDELT3}=$cdelt3;
$pdl->sethdr($h);
$pdl->wfits($out_file);

exit;






