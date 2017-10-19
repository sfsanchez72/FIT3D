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
    print "fits2ppm_hist.pl R.fits G.fits B.fits output.ppm\n";
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
print "Doing statistics...\n";
$gamma=1;


$minR=1e12;
$minG=1e12;
$minB=1e12;
$maxR=-1e12;
$maxG=-1e12;
$maxB=-1e12;
for ($i=0;$i<$nxR;$i++) {
    for ($j=0;$j<$nyR;$j++) {
	if ($minR>$a_R[$j][$i]) {
	    $minR=$a_R[$j][$i];
	}
	if ($minB>$a_B[$j][$i]) {
	    $minB=$a_B[$j][$i];
	}
	if ($minG>$a_G[$j][$i]) {
	    $minG=$a_G[$j][$i];
	}
	if ($maxR<$a_R[$j][$i]) {
	    $maxR=$a_R[$j][$i];
	}
	if ($maxB<$a_B[$j][$i]) {
	    $maxB=$a_B[$j][$i];
	}
	if ($maxG<$a_G[$j][$i]) {
	    $minG=$a_G[$j][$i];
	}
	$val_R[$i+$j*$nxR]=$a_R[$j][$i];
	$val_G[$i+$j*$nxR]=$a_G[$j][$i];
	$val_B[$i+$j*$nxR]=$a_B[$j][$i];
    }
}
$meanR=mean(@val_R);
$meanB=mean(@val_B);
$meanG=mean(@val_G);
print "R min,max=$minR,$maxR ($meanR)\n";
print "G min,max=$minG,$maxG ($meanG)\n";
print "B min,max=$minB,$maxB ($meanB)\n";
print "---------------------------\n";


$command="";
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
    print "R min,max=$minR,$maxR\n";
    print "G min,max=$minG,$maxG\n";
    print "B min,max=$minB,$maxB\n";
    $n=0;
    for ($i=0;$i<$nxR;$i++) {
	for ($j=0;$j<$nyR;$j++) {
	    $scaled_R[$i+$j*$nxR]=(255/($maxR-$minR))*($val_R[$i+$j*$nxR]-$minR);
	    $scaled_G[$i+$j*$nxG]=(255/($maxG-$minG))*($val_G[$i+$j*$nxG]-$minG);
	    $scaled_B[$i+$j*$nxB]=(255/($maxB-$minB))*($val_B[$i+$j*$nxB]-$minB);
	    $n++;
	}
    }
    pgbegin(0,"/xs",1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.2);           # Set character height
    pgsubp(3,2);
    @tr=(0,1,0,0,0,1);
    pgenv(0,$nxR,0,$nyR,1,0);
    pglabel("X","Y","Red");
    $table="ramp"; $reverse=1; $bright=0.7; $contra=0.5; $color_cont=0;
    ($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
    $nc=$pl->getdim(0);
    for ($j=0;$j<$nc;$j++) {
	$l[$j] = $pl->slice($j)->sclr;
	$r[$j] = $pr->slice($j)->sclr;
	$g[$j] = $pg->slice($j)->sclr;
	$b[$j] = $pb->slice($j)->sclr;
    }
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contra);
    pgimag(\@a_R,$nxR,$nyR,1,$nxR,1,$nyR,$minR,$maxR,\@tr);
    pgenv(0,$nxG,0,$nyG,1,0);
    pglabel("X","Y","Green");
    $table="ramp"; $reverse=1; $bright=0.7; $contra=0.5; $color_cont=0;
    ($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
    $nc=$pl->getdim(0);
    for ($j=0;$j<$nc;$j++) {
	$l[$j] = $pl->slice($j)->sclr;
	$r[$j] = $pr->slice($j)->sclr;
	$g[$j] = $pg->slice($j)->sclr;
	$b[$j] = $pb->slice($j)->sclr;
    }
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contra);
    pgimag(\@a_G,$nxG,$nyG,1,$nxG,1,$nyG,$minG,$maxG,\@tr);
    pgenv(0,$nxB,0,$nyB,1,0);
    pglabel("X","Y","Blue");
    $table="ramp"; $reverse=1; $bright=0.7; $contra=0.5; $color_cont=0;
    ($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
    $nc=$pl->getdim(0);
    for ($j=0;$j<$nc;$j++) {
	$l[$j] = $pl->slice($j)->sclr;
	$r[$j] = $pr->slice($j)->sclr;
	$g[$j] = $pg->slice($j)->sclr;
	$b[$j] = $pb->slice($j)->sclr;
    }
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contra);
    pgimag(\@a_B,$nxB,$nyB,1,$nxB,1,$nyB,$minB,$maxG,\@tr);

    pghist($n,\@val_R,$minR,$maxR,100,0);
    pghist($n,\@val_G,$minG,$maxG,100,0);
    pghist($n,\@val_B,$minB,$maxB,100,0);

    pgclos();
    pgend();

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
	system("display $output");
    }
    print "T to test the result\n";
    print "R to scales min,max in the RED\n";
    print "G to scales min,max in the GREEN\n";
    print "B to scales min,max in the BLUE\n";
    print "A to scales min,max in the All\n";
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
system("display $output");
exit;

