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
    print "USE: plot_spots.pl N.SPOT DEV SCALE\n";
    exit;
}
$ns=$ARGV[0];
$dev=$ARGV[1];
$scale=$ARGV[2];

$n=0;
$k=0;
$p=0;
open(DIR,"ls *.spots |");
while($file=<DIR>) {



    if ($file !~ "az_el") {

    chop($file);
    open(FH,"<$file");
    $line=<FH>;
    chop($line);
    @data=split(" ",$line);
    $az[$n]=$data[1]-180;
    $el[$n]=$data[2];
    $name[$n]=$file;
    $i=0;
    do {
	$line=<FH>;
	chop($line);
	@data=split(" ",$line);
	$x[$n]=$data[2];
	$y[$n]=$data[3];

	$xx[$k]=$data[2];
	$yy[$k]=$data[3];
	
	$i++;
    } while ($i<$ns);

    



    if ($az[$n]>-180) {

	$n++;
    }
    }
}
close(DIR);

$x_med=median(@x);
$y_med=median(@y);
$x_sig=sigma(@x);
$y_sig=sigma(@y);

print "$x_med+-$x_sig ; $y_med +- $y_sig\n";

for ($i=0;$i<$n;$i++) {
    $DX[$i]=($x[$i]-$x_med);
    $DY[$i]=($y[$i]-$y_med);
    $dx[$i]=$az[$i]+$scale*($x[$i]-$x_med);
    $dy[$i]=$el[$i]+$scale*($y[$i]-$y_med);
#    print "$i $dx[$i] $dy[$i] $x[$i] $y[$i]\n";
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv(-180,180,0,90.1,0,0);
pglabel("Azimuth","Elevation","");
pgsci(3);
pgline($n,\@az,\@el);
pgsci(1);
pgsch(0.7);
for ($i=0;$i<$n;$i++) {
    pgarro($az[$i],$el[$i],$dx[$i],$dy[$i]);
}
pgsch(1.7);
pgsci(2);
pgpoint($n,\@az,\@el,1);

#pgline($n,\@dx,\@dy,16);
pgclose;
pgend;



$P[0]=0;
#print "$name[0] 0 $p\n";
for ($i=1;$i<$n;$i++) {
    if (abs($el[$i]-$el[$i-1])>5) {
	$p++;
    }
    $P[$i]=$p;
#    print "$name[$i] $i $p\n";
}

if ($dev !~ "PS") {
    print "Press Enter\n"; <stdin>;
}

pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv(-181,181,-1.5,1.5,0,0);
pglabel("Azimuth","Delta X,Y","");
pgsch(1.7);
pgsci(2);
pgpoint($n,\@az,\@DX,16);
pgsci(5);
pgpoint($n,\@az,\@DY,22);
pgclose;
pgend;

if ($dev !~ "PS") {
    print "Press Enter\n"; <stdin>;
}



pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv(-185,185,-1.5,1.5,0,0);
pglabel("Azimuth","Delta X,Y","");
pgsch(1.7);

open(OUT,">az_el.spots");
for ($i=0;$i<$n;$i++) {
    $color=$P[$i]+1;
    $XX=$az[$i];
    $YY=$DX[$i];
 #   print "$XX $YY ";
    pgsci($color);
    pgpoint(1,[$XX],[$YY],16);
    $YY=$DY[$i];
    pgpoint(1,[$XX],[$YY],22);
  #  print "$YY \n";

    $DD=sqrt($DX[$i]**2+$DY[$i]**2);
    print OUT "$az[$i] $el[$i] $DD\n";
    print "$name[$i] $az[$i] $el[$i] $DX[$i] $DY[$i] $DD\n";
}
close(OUT);
$x_sig=apr($x_sig);
$y_sig=apr($y_sig);
pgsch(1.4);
pgsci(1);
pgptxt(-172,1.15,0,0,"\\gs(\\gDx)=$x_sig pix");
pgptxt(-172,0.95,0,0,"\\gs(\\gDy)=$y_sig pix");

pgclose;
pgend;

exit;

