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
    print "USE: res2var.pl INPUT_RES.FITS OUTPUT_VAR.fits radius\n";
    exit;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$pdl=rfits($infile);
($nx,$ny,$nz)=$pdl->dims();
$mpdl=zeroes($nx,$ny,$nz);
$pdl->inplace->setnantobad;
$pdl->inplace->setbadtoval(0);
$r=$ARGV[2];



$cdelt=$pdl->hdr->{"CDELT"};
if ($cdelt<1) {
    $cdelt=1;
}
$box=$cdelt*5;
for ($i=0;$i<$nx;$i++) {
    my $t=$pdl->slice("($i),,");
    @n=$t->dims;
    my $soft=$mpdl->slice("($i),,");
    $t->inplace->abs;
    $soft .= med2df($t,1,$box,{Boundary=>Reflect});
    print "$i/$nx\n";
}


for ($i=0;$i<$nz;$i++) {
    my $t=$pdl->slice(",,($i)");
    @n=$t->dims;
    my $soft=$mpdl->slice(",,($i)");
    $soft .= med2df($soft,$r,$r,{Boundary=>Reflect});
    $soft .= $soft;
    $soft .= $soft->setnantobad;
    $soft .= $soft->setbadtoval(10e20);
}

$mpdl=$mpdl/$r;
$h=$pdl->gethdr;
$mpdl->sethdr($h);
$mpdl->wfits($outfile);
$a=sqrt(abs($mpdl));
$a->wfits("sigma.fits");

exit;

