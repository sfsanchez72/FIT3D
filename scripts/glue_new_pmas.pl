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


if ($#ARGV<0) {
    print "USE: glue_new_pmas.pl PREFIX\n";
    exit;
}

$prefix=$ARGV[0];
$outfile=$prefix.".fits";
$call="rm -f ".$outline;
system($call);
$in_a=$prefix."a.fits";
$in_b=$prefix."b.fits";
$in_c=$prefix."c.fits";
$in_d=$prefix."d.fits";

$pdl_a=rfits($in_a);
$pdl_b=rfits($in_b);
$pdl_c=rfits($in_c);
$pdl_d=rfits($in_d);

$hdr=$pdl_a->gethdr;


$gain_a=$pdl_a->hdr->{GAIN};
$gain_b=$pdl_b->hdr->{GAIN};
$gain_c=$pdl_c->hdr->{GAIN};
$gain_d=$pdl_d->hdr->{GAIN};

$pdl_a=$pdl_a*$gain_a;
$pdl_b=$pdl_b*$gain_b;
$pdl_c=$pdl_c*$gain_c;
$pdl_d=$pdl_d*$gain_d;


$pdl_out=zeroes(float,2048,2056);

#30:1030,10:30
($mean_a,$rms_a,$median_a,$min_a,$max_a) = stats($pdl_a->slice('15:23,15:1022'));
($mean_b,$rms_b,$median_b,$min_b,$max_b) = stats($pdl_b->slice('15:23,15:1022'));
($mean_c,$rms_c,$median_c,$min_c,$max_c) = stats($pdl_c->slice('15:23,15:1022'));
($mean_d,$rms_d,$median_d,$min_d,$max_d) = stats($pdl_d->slice('15:23,15:1022'));

$s_pdl_a=$pdl_a->slice('25:1048,0:1027');#-$median_a;
$s_pdl_b=$pdl_b->slice('25:1048,0:1027');#-$median_b;
$s_pdl_c=$pdl_c->slice('25:1048,0:1027');#-$median_c;
$s_pdl_d=$pdl_d->slice('25:1048,0:1027');#-$median_d;


#$pdl_b=$pdl_a->slice('-*:*');

($nx,$ny)=$s_pdl_a->dims();

#A
$nx0=0;
$nx1=$nx-1;
$ny0=2055;
$ny1=2055-$ny+1;
$sec="".$nx0.":".$nx1.",".$ny0.":".$ny1;
$t=$pdl_out->slice($sec);
$t.=$s_pdl_a-$median_a;

$nx0=$nx;
$nx1=$nx+$nx-1;
$ny0=2055;
$ny1=2055-$ny+1;
$sec="".$nx0.":".$nx1.",".$ny0.":".$ny1;
$t=$pdl_out->slice($sec);
$nx1=$nx-1;
$t.=$s_pdl_b->slice("$nx1:0,")-$median_b;


$nx0=0;
$nx1=$nx-1;
$ny0=2055-1028;
$ny1=2055-1028-$ny+1;
$sec="".$nx0.":".$nx1.",".$ny0.":".$ny1;
$t=$pdl_out->slice($sec);
$ny1=$ny-1;
$t.=$s_pdl_d->slice(",$ny1:0")-$median_d;

$nx0=$nx;
$nx1=$nx+$nx-1;
$ny0=2055-1028;
$ny1=2055-1028-$ny+1;
$sec="".$nx0.":".$nx1.",".$ny0.":".$ny1;
$t=$pdl_out->slice($sec);
$ny1=$ny-1;
$nx1=$nx-1;
$t.=$s_pdl_c->slice("$nx1:0,$ny1:0")-$median_c;


$pdl_out->sethdr($hdr);
$pdl_out->hdr->{BITPIX}=-32;
$pdl_out->hdr->{NAXIS1}=2048;
$pdl_out->hdr->{NAXIS2}=2056;
$pdl_out->wfits($outfile);


exit;


#C
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$s_pdl_c->at($nx-1-$i,$ny-1-$j)-$median_c;
	set($pdl_out,$nx+$i,2055-(1028+$j),$val);
    }
}





