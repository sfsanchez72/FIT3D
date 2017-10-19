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



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<1) {
    print "USE: check_trace.pl TRACE.fits y_width\n";
    exit;
}

$trace_file=$ARGV[0];
$y_width=$ARGV[1];

print "Reading file $trace_file ";
$nax=read_naxes($trace_file);   
@naxis=@$nax;
@peak_y_max_new=read_img($trace_file);
print "Done\n";
$nx=$naxis[0];
$ny=$naxis[1];


pgbegin(0,"/xs",1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.2);           # Set character height
for ($k=0;$k<$ny;$k++) {
    pgenv(0,$nx,$peak_y_max_new[$k][int($nx/2)]-$y_width,$peak_y_max_new[$k][int($nx/2)]+$y_width,0,0);
    pglabel("X-axis","Y-Axis","$k/$ny");
    for ($i=0;$i<$nx;$i++) {
	pgsci(4);
	$x=$i;
	$y=$peak_y_max_new[$k][$i];
	pgpoint(1,[$x],[$y],5);
    }
    pgsci(1);
    pgclose;
#    print "Press Enter"; <stdin>; 
}
pgend;

exit;
