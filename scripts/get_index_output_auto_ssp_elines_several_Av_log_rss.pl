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

#use POSIX;



$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<3) {
    print "USE:get_index_rss_montecarlo.pl output_auto_ssp_elines_several.RSS.cube.fits NSIM REDSHIFT DEV\n";    
    exit;
}

$spec_file=$ARGV[0];
$NSIM=$ARGV[1];
$z=$ARGV[2];
$dev=$ARGV[3];

$pdl=rfits("$spec_file");
($nx,$ny,$nz)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL3"};
$cdelt=$pdl->hdr->{"CDELT3"};
$crpix=$pdl->hdr->{"CRPIX3"};


for ($j=0;$j<$ny;$j++) {
    open(FH,">get_index_rss.spec");
    open(FH1,">get_index_noise.spec");
    for ($i=0;$i<$nx;$i++) {
	$wave=$crval+$cdelt*($i+1-$crpix);
	$flux_org=$pdl->at($i,$j,0);
	$flux_res=$pdl->at($i,$j,4);
#	print "$wave $flux_org $flux_res $flux_res1 $flux_res2\n";
	if ($flux_res==0) {
	    $flux_res=sqrt(abs($flux_org));
	}
	print FH "$i $wave $flux_org\n";
	print FH1 "$i $wave $flux_res\n";
    }
    close(FH);
    close(FH1);
    $call="get_index_montecarlo.pl get_index_rss.spec get_index_noise.spec ".$NSIM." ".$z." ".$dev;
    system($call);
}
exit;
