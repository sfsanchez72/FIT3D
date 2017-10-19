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

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";
$ENV{'PGPLOT_ENVOPT'}="V";


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";
if ($#ARGV<2) {
    print "USE: flux_filters.pl LIST_OF_FILTERS spectrum redshift \n";
    exit;
}

$C=299792.458;

$filters=$ARGV[0];
$input=$ARGV[1];
$redshift=$ARGV[2];
$vel_mod=$redshift*$C;


open(FH,"<$filters");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$filter=$data[0];
	$call="flux_filter.pl ".$filter." ".$input." ".$redshift." > junk.junk";
	system($call);
#	print "$call\n";
#	exit;
	open(CALL,"<flux_filter.out");

	while ($out=<CALL>) {
	    chop($out);
	    if ($out!~ "some") {
		print "$out\n";
	    }
	}
	close(CALL);
    }
}
close(FH);

exit;



