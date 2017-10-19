#!/usr/bin/perl
#
#
#
use PDL;
use PDL::Core;
use PDL::Graphics::LUT;


use Astro::FITS::CFITSIO;
use PGPLOT;
use Carp;

use PDL::IO::Pic;



$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";



if ($#ARGV<5) {
    print "fits2ppm.pl R.fits G.fits B.fits output.ppm min max\n";
    exit;
}

$Rfile=$ARGV[0];
$Gfile=$ARGV[1];
$Bfile=$ARGV[2];
#$max=$ARGV[3];
$output=$ARGV[3];
$min=$ARGV[4];
$max=$ARGV[5];

print "Start Reading ";
$pdl_R=rfits($Rfile);
$hdr=$pdl_R->hdr;
#$hdr=rfits($Rfile,{data=>0});
$nxR = $hdr->{NAXIS1};
$nyR = $hdr->{NAXIS2};
print "DIMS=[$nxR,$nyR]\n";
$stack = zeroes(3,$nxR,$nyR);
$pdl_R=rfits($Rfile);
$pdl_G=rfits($Gfile);
$pdl_B=rfits($Bfile);
($nxR,$nyR) = $pdl_R->dims;
($nxG,$nyG) = $pdl_G->dims;
($nxB,$nyB) = $pdl_B->dims;
#($min,$max) = $pdl_R->minmax;

print "Done\n";
#print "Doing statistics...\n";
$gamma=1;
	$t = $stack->slice('(0),:,:');
	$t .=$R;

	$t = $stack->slice('(1),:,:');
	$t .=$G;
	
	$t = $stack->slice('(2),:,:');
	$t .=$B;

print "MINMAX=$min $max\n";	
$R=byte(zeroes($nxR,$nyR));
$G=byte(zeroes($nxG,$nyG));
$B=byte(zeroes($nxB,$nyB));
$R=byte((255/($max-$min))*($pdl_R-$min));
$G=byte((255/($max-$min))*($pdl_G-$min));
$B=byte((255/($max-$min))*($pdl_B-$min));



$t = $stack->slice('(0),:,:');
$t .=$R;

$t = $stack->slice('(1),:,:');
$t .=$G;

$t = $stack->slice('(2),:,:');
$t .=$B;

$stack->wpic($output);

exit;

