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

print "Reading $Rfile\n";
my $naxes=read_naxes($Rfile);
($nxR,$nyR) = @$naxes;
@a_R=read_img($Rfile);
print "Reading $Gfile\n";
my $naxes=read_naxes($Gfile);
($nxG,$nyG) = @$naxes;
@a_G=read_img($Gfile);
print "Reading $Bfile\n";
my $naxes=read_naxes($Bfile);
($nxB,$nyB) = @$naxes;
@a_B=read_img($Bfile);
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
	open(FH,">$output");
	print FH "P3\n";
	print FH "# $output\n";
	print FH "$nyR $nxR\n";
	print "$nxR $nyR\n";
	print FH "256\n";
	for ($i=0;$i<$nxR;$i++) {
	    for ($j=0;$j<$nyR;$j++) {
		$r=int($scaled_R[$i+$j*$nxR]);
		$g=int($scaled_G[$i+$j*$nxG]);
		$b=int($scaled_B[$i+$j*$nxB]);	
		if ($r<0) {
		    $r=0;
		}
		if ($g<0) {
		    $g=0;
		}
		if ($b<0) {
		    $b=0;
		}
		if ($r>255) {
		    $r=255;
		}
		if ($g>255) {
		    $g=255;
		}
		if ($b>255) {
		    $b=255;
		}

		print FH "$r $g $b\n";
	    }
	}
	close(FH);
	system("xv $output");
    }

    print "R min,max=$minR,$maxR\n";
    print "G min,max=$minG,$maxG\n";
    print "B min,max=$minB,$maxB\n";
    if ($command ne "T") {
	$n=0;
	for ($i=0;$i<$nxR;$i++) {
	    for ($j=0;$j<$nyR;$j++) {
		$scaled_R[$i+$j*$nxR]=(255/($maxR-$minR))*($a_R[$j][$i]-$minR);
		$scaled_G[$i+$j*$nxG]=(255/($maxG-$minG))*($a_G[$j][$i]-$minG);
		$scaled_B[$i+$j*$nxB]=(255/($maxB-$minB))*($a_B[$j][$i]-$minB);
#	    print "$val_R[$i+$j*nxR] $minR $maxR $scaled_R[$i+$j*$nxR]\n";
	    $n++;
	    }
	}
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




open(FH,">$output");
print FH "P3\n";
print FH "# $output\n";
print FH "$nyR $nxR\n";
print "$nxR $nyR\n";
print FH "256\n";
for ($i=0;$i<$nxR;$i++) {
    for ($j=0;$j<$nyR;$j++) {
	$r=int($scaled_R[$i+$j*$nxR]);
	$g=int($scaled_G[$i+$j*$nxG]);
	$b=int($scaled_B[$i+$j*$nxB]);	
	if ($r<0) {
	    $r=0;
	}
	if ($g<0) {
	    $g=0;
	}
	if ($b<0) {
	    $b=0;
	}
	if ($r>255) {
	    $r=255;
	}
	if ($g>255) {
	    $g=255;
	}
	if ($b>255) {
	    $b=255;
	}

	print FH "$r $g $b\n";
    }
}
close(FH);
system("xv $output");
exit;

