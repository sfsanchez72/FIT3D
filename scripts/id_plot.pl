#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
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

$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

#$ENV{"PGPLOT_FONT"}="/work/sanchez/sda2/pgplot/grfont.dat";


if ($#ARGV<2) {
    print "USE: id_plot.pl INPUT_FILE DEVICE emission_lines.txt [MIN MAX] [WMIN WMAX]\n";
    exit;
}

$input=$ARGV[0];
$dev=$ARGV[1];
$efile=$ARGV[2];
$def=0;
if ($#ARGV==4) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $def=1;
}

$ne=0;
open(FH,"<$efile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$we[$ne]=$data[0];
	$ide[$ne]=$data[1]." ".$data[0];
	$ne++;
    }
}
close(FH);
print "$ne emission lines found\n";

$y_min=1e12;
$y_max=-1e12;
$n=0;
open(FH,"<$input");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $flux[$n]=$data[2];#*(8281/0.65)*5.6410*1e-8;
    if ($flux[$n]>0) {
	$lflux[$n]=log10($flux[$n]);
    } else {
	$lflux[$n]=-1000;
    }
    $wave[$n]=$data[1];
    if ($flux[$n]>$y_max) {
	$y_max=$flux[$n];
    }
    if ($flux[$n]<$y_min) {
	$y_min=$flux[$n];
    }
    $n++;
}
close(FH);
$wmin=$wave[0];
$wmax=$wave[$n-1];
if ($#ARGV==6) {
    $min=$ARGV[3];
    $max=$ARGV[4];
    $wmin=$ARGV[5];
    $wmax=$ARGV[6];
    $def=1;
}

print "$#ARGV $min $max $wmin $wmax\n";
if ($def==1) {
    $y_min=$min;    
    $y_max=$max;
}

pgbegin(0,$dev,1,1);
#pgpap(14.0,0.4);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv($wmin,$wmax,$y_min,$y_max,0,0);
#pgenv($wmin,$wmax,$y_min,$y_max,0,20);
pglabel("Wavelength (\\A)","Flux ","");
pgsls(1);
#pgline($n,\@wave,\@lflux);
pgline($n,\@wave,\@flux);
pgsch(1.1);
for ($i=0;$i<$ne;$i++) {
#    pgptxt($we[$i],0.8*$y_max,90,0,$ide[$i]);
    pgsci(8);

    if ($i==(2*int($i/2))) {
	$pos=0.7*$y_max;
    } else {
	$pos=0.6*$y_max;
    }
#    pgptxt($we[$i],0.5*$y_max,0,0,"$ide[$i]");
    pgptxt($we[$i]+3,$pos,90,0,"$ide[$i]");
    
    pgsci(14);
    
    pgline(2,[$we[$i],$we[$i]],[0.9*$pos,0.7*$pos]);
#    print "$we[$i] $pos $ide[$i]\n";
}
pgclose;
pgend;

exit;
