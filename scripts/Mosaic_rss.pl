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

if ($#ARGV<2) {
    print "USE: Mosaic_rss.pl CONFIG.txt RSS.FITS POSITION_TABLE.txt [NO_FITS]\n";
    print "CONFIG: infile in_position_table dx dy\n";
    exit;
}

$config=$ARGV[0];
$outfile=$ARGV[1];
$outpt=$ARGV[2];
$no_fits=0;
if ($#ARGV==3) {
    $no_fits=1;
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
    $inpt[$nf]=$data[1];
    $dx[$nf]=$data[2];
    $dy[$nf]=$data[3];
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

#
# Reading Position tables:
#
open(OUT,">$outpt");
$n=0;
for ($j=0;$j<$nf;$j++) {
    $pt_file=$inpt[$j];
    print "Adding $pt_file, ";
    open(FH,"<$pt_file");
    if ($j==0) {
	$line=<FH>;
	print OUT "$line";
    } else {
	$line=<FH>;
    }
    while($line=<FH>) {
	chop($line);
        @data=split(" ",$line);
	$id=$data[0]+$j*$n;
	if ($j==0) {
	    $n++;
	}
	$x=$data[1]+$dx[$j];
	$y=$data[2]+$dy[$j];
	$type=$data[3];
	print OUT "$id $x $y $type\n";
    }
    close(FH);
    print "done\n";
    if ($no_fits==0) {
	$a_in=rfits($infile[$j]);
	if ($j==0) {
	    $h=$a_in->gethdr;
	    ($nx,$ny)=$a_in->dims;
	    $ny_out=$nf*$ny;
	    $a_out=zeroes($nx,$ny_out);
	}
	for ($jj=0;$jj<$ny;$jj++) {
	    $jjj=$jj+$j*$ny;
	    for ($i=0;$i<$nx;$i++) {	    
		$val=$a_in->at($i,$jj);
		set($a_out,$i,$jjj,$val);
	    }
	}
    }
    print "file $j added ($nx,$ny_out)\n";
}
close(OUT);
if ($no_fits==0) {
    $h->{NAXIS2}=$ny_out;
    $a_out->sethdr($h);
    $a_out->wfits($outfile);
}
print "DONE\n";
exit;


