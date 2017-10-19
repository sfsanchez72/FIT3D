#!/usr/bin/perl
#
#
#


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");


if ($#ARGV<6) {
    print "fits2ppm.pl R.fits G.fits B.fits output.ppm scale_R scale_G scale_B\n";
    exit;
}

$Rfile=$ARGV[0];
$Gfile=$ARGV[1];
$Bfile=$ARGV[2];
#$max=$ARGV[3];
$output=$ARGV[3];
$scR=$ARGV[4];
$scG=$ARGV[5];
$scB=$ARGV[6];

print "Start Reading\n";
my $naxes=read_naxes($Rfile);
($nxR,$nyR) = @$naxes;
@a_R=read_img($Rfile);
my $naxes=read_naxes($Gfile);
($nxG,$nyG) = @$naxes;
@a_G=read_img($Gfile);
my $naxes=read_naxes($Bfile);
($nxB,$nyB) = @$naxes;
@a_B=read_img($Bfile);
print "Start Reading\n";

open(FH,">$output");
print FH "P3\n";
print FH "# $output\n";
print FH "$nyR $nxR\n";
print "$nxR $nyR\n";
print FH "256\n";
for ($i=0;$i<$nxR;$i++) {
    for ($j=0;$j<$nyR;$j++) {
	$r=int($a_R[$j][$i]*$scR);
	$g=int($a_G[$j][$i]*$scG);
	$b=int($a_B[$j][$i]*$scB);
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
exit;

