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





$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";



if ($#ARGV<3) {
    print "fits2ppm.pl R.fits G.fits B.fits output.ppm\n";
    exit;
}

$Rfile=$ARGV[0];
$Gfile=$ARGV[1];
$Bfile=$ARGV[2];
#$max=$ARGV[3];
$output=$ARGV[3];

print "Start Reading ";
$hdr=rfits($Rfile,{data=>0});
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
($min,$max) = $pdl_R->minmax;

print "Done\n";
#print "Doing statistics...\n";
$gamma=1;


$command="A";
while ($command ne "Q") {
    if ($command eq "A") {
	print "minR maxR minG maxG minB maxB=";
	$minmax=<STDIN>;
	chop($minmax);
	($minR,$maxR,$minG,$maxG,$minB,$maxB)=split(" ",$minmax);
    }
    if ($command eq "R") {
	print "min max=";
	$minmax=<STDIN>;
	chop($minmax);
	($minR,$maxR)=split(" ",$minmax);
    }
    if ($command eq "G") {
	print "min max=";
	$minmax=<STDIN>;
	chop($minmax);
	($minG,$maxG)=split(" ",$minmax);
    }
    if ($command eq "B") {
	print "min max=";
	$minmax=<STDIN>;
	chop($minmax);
	($minB,$maxB)=split(" ",$minmax);
    }

    if ($command eq "T") {
	$t = $stack->slice('(0),:,:');
	$t .=$R;

	$t = $stack->slice('(1),:,:');
	$t .=$G;
	
	$t = $stack->slice('(2),:,:');
	$t .=$B;
	
	$stack->wpic($output);

	system("xv $output");
    }

    print "R min,max=$minR,$maxR\n";
    print "G min,max=$minG,$maxG\n";
    print "B min,max=$minB,$maxB\n";
    if ($command ne "T") {
	$R=byte(zeroes($nxR,$nyR));
	$G=byte(zeroes($nxG,$nyG));
	$B=byte(zeroes($nxB,$nyB));
	$R=byte((255/($maxR-$minR))*($pdl_R-$minR));
	$G=byte((255/($maxG-$minG))*($pdl_G-$minG));
	$B=byte((255/($maxB-$minB))*($pdl_B-$minB));

	print "$R\n";
    }

    print "R to scales min,max in the RED\n";
    print "G to scales min,max in the GREEN\n";
    print "B to scales min,max in the BLUE\n";
    print "A to scales min,max in the All\n";
    print "T to test the result\n";
    print "Q to quit\n";
    $command=<STDIN>;
    chop($command);
}



	$t = $stack->slice('(0),:,:');
	$t .=$R;

	$t = $stack->slice('(1),:,:');
	$t .=$G;
	
	$t = $stack->slice('(2),:,:');
	$t .=$B;
	
	$stack->wpic($output);

exit;

