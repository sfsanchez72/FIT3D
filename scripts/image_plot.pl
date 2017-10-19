#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
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

use PDL::Core;
use PDL::Graphics::LUT;
use Carp;

$ENV{'PGPLOT_ENVOPT'}="IV";


#$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";
#$galfit="/home/sanchez/sda1/galfit/galfit";
#$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<1) {
    print "USE: image_plot.pl INPUT_FILE.FITS DEVICE [min max]\n";
    exit;
}

$input=$ARGV[0];
$dev=$ARGV[1];

$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;
#($min,$max)=minmax($pdl);
#$table="ramp"; $reverse=0; 
$table="ramp"; $reverse=1; 
#$table="heat"; $reverse=0; 
@map=list($pdl);
$median=median($pdl);
$new_pdl=($pdl-$median)**2;
$sigma=median($new_pdl);
#$sigma=sigma(@map);
#
$min=$median-$sigma;
$max=$median+3*$sigma;

if ($#ARGV==3) {
    $min=$ARGV[2];
    $max=$ARGV[3];
}
print "$min $max $median $sigma\n";

@tr=(0,1,0,0,0,1);

pgbeg(0,$dev,1,1);
pgscf(2.0);
pgsch(0.8);
($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
$nc=$pl->getdim(0);
for ($j=0;$j<$nc;$j++) {
    $l[$j] = $pl->slice($j)->sclr;
    $r[$j] = $pr->slice($j)->sclr;
    $g[$j] = $pg->slice($j)->sclr;
    $b[$j] = $pb->slice($j)->sclr;
}
$bright=0.9;
$contrast=0.7;
if ($#ARGV==5) {
    $min=$ARGV[2];
    $max=$ARGV[3];
    $bright=$ARGV[4];
    $contrast=$ARGV[5];
}
pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);
pgenv(0,$nx,0,$ny,1,0);
#pgpap(6.0,1.0);

#pgsch(0.7);           # Set character height
#pgsch(1.2);

pglabel("X-pix","Y-pix","");
pgimag(\@map,$nx,$ny,1,$nx,1,$ny,$min,$max,\@tr);
pgsci(1);
#pgwedg("RI",0.0,4.5,$min,$max,"");
#pgbox("SBC",0,0,"SBC",0,0);


pgclos();
pgend();



exit;
