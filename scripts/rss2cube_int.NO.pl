#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module
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

use PDL::NiceSlice;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<6) {
    print "USE: rss2cube_int.pl INPUT_RSS position_table.txt dx OUTPUT.cube.fits alpha BOX SCALE\n";
    print "type:\n";
    print "weight = exp(-0.5*(dist/(SCALE*dx))**alpha)\n";
    print "Recommended alpha=2, BOX=5, SCALE = 1\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$dx=$ARGV[2];
$output_rss=$ARGV[3];
$alpha=$ARGV[4];
$A=$ARGV[5];
$SCALE=$ARGV[6];



$n=0;
open(PT,"<$pos_table");
$line=<PT>;
chop($line);
@data=split(" ",$line);
$R=$data[1];
$nx_min=1e12;
$ny_min=1e12;
$nx_max=-1e12;
$ny_max=-1e12;

$x_min=1e12;
$y_min=1e12;
$x_max=-1e12;
$y_max=-1e12;


while($line=<PT>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $xx[$n]=$data[1];
    $yy[$n]=$data[2];
    if ($xx[$n]>$x_max) {
	$x_max=$xx[$n];
    }
    if ($xx[$n]<$x_min) {
	$x_min=$xx[$n];
    }
    if ($yy[$n]>$y_max) {
	$y_max=$yy[$n];
    }
    if ($yy[$n]<$y_min) {
	$y_min=$yy[$n];
    }
    $n++;
}
close(PT);
$x_min=$x_min-2*$R;
$x_max=$x_max+2*$R;
$y_min=$y_min-2*$R;
$y_max=$y_max+2*$R;
$nxc=int(($x_max-$x_min)/$dx);
$nyc=int(($y_max-$y_min)/$dx);


$input=rfits($input_rss);
$nx=$input->hdr->{NAXIS1};
$ny=$input->hdr->{NAXIS2};
#($nx,$ny)=$input->dims();
$nzc=$nx;

#print "($nxc,$nyc,$nzc) $nx $ny\n";
$output=zeroes($nxc,$nyc,$nzc);
$mask=zeroes($nxc,$nyc);
$pdl_w=zeroes($nxc,$nyc);
$input_mask=zeroes($nx,$ny);
$output_mask=zeroes($nxc,$nyc,$nzc);

for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$input->at($i,$j);
	if ($val!=0) {
	    set($input_mask,$i,$j,1);
	} else {
	    set($input_mask,$i,$j,0);
	}
    }
}



for ($ii=0;$ii<$n;$ii++) {
    $ic=int(($xx[$ii]-$x_min)/$dx);
    $jc=int(($yy[$ii]-$y_min)/$dx);
    $ic_val=(($xx[$ii]-$x_min)/$dx);
    $jc_val=(($yy[$ii]-$y_min)/$dx);
    $spec=$input->slice(",($ii)");
    $spec_mask=$input_mask->slice(",($ii)");

    $D=$SCALE*$dx;#sqrt($nxc**2+$nyc*+2);


#
# Influence area
#    
    $i_min=$ic-$A*$dx;
    $i_max=$ic+$A*$dx;
    $j_min=$jc-$A*$dx;
    $j_max=$jc+$A*$dx;

    if ($i_min<0) {
	$i_min=0;
    }
    if ($j_min<0) {
	$j_min=0;
    }
    if ($i_max>=$nxc) {
	$i_max=$nxc-1;
    }
    if ($j_max>=$nyc) {
	$j_max=$nyc-1;
    }

    for ($i=$i_min;$i<$i_max;$i++) {
	for ($j=$j_min;$j<$j_max;$j++) {
	    $dist=sqrt(($i-$ic_val)**2+($j-$jc_val)**2);
	    $w=exp(-0.5*($dist/$D)**$alpha);
	    $w_in=$output_mask->slice("($i),($j),");
	    my $spec_mask_now=$spec_mask*$w;
	    $w_in .= $w_in+$spec_mask_now;
	    my $t = $output->slice("($i),($j),");	    
	    $a=$spec*$spec_mask_now;
	    $t .= $t+$a;
	    if ($dist<2*$D) {
		set($mask,$i,$j,1);
	    }
	}
    }
#    print "$ii/$n\n";
}

$output=$output/$output_mask;
$output=$output*$mask;

#
# Final mask
#




$crval=$input->hdr->{CRVAL1};
$cdelt=$input->hdr->{CDELT1};
$output->hdr->{CRPIX1}=1;
$output->hdr->{CRVAL1}=$x_min;
$output->hdr->{CDELT1}=$dx;
$output->hdr->{CRPIX2}=1;
$output->hdr->{CRVAL2}=$y_min;
$output->hdr->{CDELT2}=$dx;
$output->hdr->{CRPIX3}=1;
$output->hdr->{CRVAL3}=$crval;
$output->hdr->{CDELT3}=$cdelt;
$output->wfits($output_rss);


exit;

