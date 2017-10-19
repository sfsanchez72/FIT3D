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



#$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<6) {
    print "USE: fiber_trans_cube_.pl SKYFLAT.fits FIBERFLAT.fits NX_min_E NX_max_E NX_min_con NX_max_con N.BOX\n";
    exit;
}

$infile=$ARGV[0];
$out_file=$ARGV[1];
$nx_min=$ARGV[2];
$nx_max=$ARGV[3];
$nxc_min=$ARGV[4];
$nxc_max=$ARGV[5];
$nbox=$ARGV[6];


#$fiber_flat=$ARGV[2];


$pdl_in=rfits($infile);
($nx,$ny,$nz)=$pdl_in->dims();
$h=$pdl_in->gethdr();

$pdl_out=ones($nx,$ny,$nz);

$pdl_tmp=zeroes(($nx/$nbox),($ny/$nbox));
$ii=0;
$jj=0;
$sumA=0;
$nA=0;
for ($i=0;$i<$nx;$i=$i+$nbox) {
    for ($j=0;$j<$ny;$j=$j+$nbox) {
	$ie=$i+$nbox-1;
	$je=$j+$nbox-1;
#	print "$i,$ie $j,$je $nx,$ny\n";
	my $pdl_eline=$pdl_in->slice("$i:$ie,$j:$je,$nx_min:$nx_max");
	my $pdl_cont=$pdl_in->slice("$i:$ie,$j:$je,$nxc_min:$nxc_max");
	$pdl_eline=$pdl_eline;#-$pdl_cont;
	$sum = median($pdl_eline);
#	print "$sum\n";
	$sumA=$sumA+$sum;
	$nA++;
	for ($ii=$i;$ii<$ie+1;$ii++) {
	    for ($jj=$j;$jj<$je+1;$jj++) {		
		my $t=$pdl_out->slice("($ii),($jj),:");
		$t .= $t*$sum;
	    }
	}
    }
}
$sumA=$sumA/$nA;

$pdl_out=$pdl_out/$sumA;
$pdl_out->sethdr($h);
$pdl_out->wfits($out_file);

exit;
