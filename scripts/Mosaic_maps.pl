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

if ($#ARGV<1) {
    print "USE: Mosaic_maps.pl CONFIG.txt OUTPUT [THRESHOLD]\n";
    exit;
}

$config=$ARGV[0];
$outfile=$ARGV[1];
$limit=-1e40;
if ($#ARGV==2) {
    $limit=$ARGV[2];
}
$nf=0;
$dxmin=1e12;
$dxmax=-1e12;
$dymin=1e12;
$dymax=-1e12;
open(FH,"<$config");
while($line=<FH>) {
#    chop($line);
    @data=split(" ",$line);
    $infile[$nf]=$data[0];
    $dx[$nf]=$data[1];
    $dy[$nf]=$data[2];
    if ($dxmin>$dx[$nf]) {
	$dxmin=$dx[$nf];
    }
    if ($dxmax<$dx[$nf]) {
	$dxmax=$dx[$nf];
    }
    if ($dymin>$dy[$nf]) {
	$dymin=$dy[$nf];
    }
    if ($dymax<$dy[$nf]) {
	$dymax=$dy[$nf];
    }
    $nf++;
}
close(FH);
print "$nf files\n";
$nax=read_naxes($infile[0]);   
@naxis=@$nax;
$nx=$naxis[0];
$ny=$naxis[1];
#$nz=$naxis[2];
print "[$nx,$ny,$nz]\n";
@tr=read_img_headers($infile[0],["CRPIX1","CRVAL1","CDELT1","CRPIX2","CRVAL2","CDELT2"]);


$xmin=int($dxmin);
$xmax=int($nx+$dxmax);
$ymin=int($dymin);
$ymax=int($ny+$dymax);

$nx_new=$xmax-$xmin;
$ny_new=$ymax-$ymin;

print "$xmin $xmax $ymin $ymax $nx_new $ny_new\n";

for ($j=0;$j<$ny_new;$j++) {
    for ($i=0;$i<$nx_new;$i++) {
	$out_array[$j][$i]=0;
    }
}


for ($f=0;$f<$nf;$f++) {
    print "Reading file $infile[$f], ";  
    my $pdl_array=rfits($infile[$f]);
    $nx_start=$dx[$f]-$dxmin;
    $ny_start=$dy[$f]-$dymin;
    print "added at $nx_start,$ny_start\n";
    for ($j=0;$j<$ny;$j++) {
	for ($i=0;$i<$nx;$i++) {
	    $val=$pdl_array->at($i,$j);
	    if ($val>$limit) {
		$out_array[$j+$ny_start][$i+$nx_start]+=$val;
		$sum[$j+$ny_start][$i+$nx_start]+=1;
	    }
	}
    }
}

for ($j=0;$j<$ny_new;$j++) {
    for ($i=0;$i<$nx_new;$i++) {
	if ($sum[$j][$i]>0) {
	    $out_array[$j][$i]/=$sum[$k][$j][$i];
	}
    }
}

$pdl_out=pdl(@out_array);
print "Writting $outfile\n";
$pdl_out->wfits($outfile);
print "DONE\n";
exit;


