#!/usr/bin/perl
#
# This programs creates a set of fits files
# extracting the information from the line-fitting results.
#


use Statistics::OLS;
use Math::Stat;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);

use Math::Approx;


use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Core;
use PDL::Graphics::LUT;
use Carp;
$ENV{PGPLOT_FOREGROUND} = "black";
$ENV{PGPLOT_BACKGROUND} = "white";





if ($#ARGV<4) {
    print "USE: plot_slice.pl in_slice.txt in_slice2.txt DEVICE SHAPE SIZE [MARK] [MIN MAX] [BRIGHT CONTRAST LUT_TABLE INVERSE(0/1)]\n";
    exit;
}

$slice=$ARGV[0];
$slice2=$ARGV[1];
$dev=$ARGV[2];
$shape=$ARGV[3];
$size=$ARGV[4];

$nmark=-1;
if ($#ARGV>=5) {
    $nmark=$ARGV[5];
}
#$cut=$ARGV[1];
#$outfile=$ARGV[2];


$ns=0;
$x_min=10e10;
$x_max=-10e10;
$y_min=10e10;
$y_max=-10e10;
$flux_max=-1e20;
$flux_min=1e20;
$dpix=10e10;
open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	$flux[$ns]=$data[3];

	if ($x_min>$x[$ns]) {
	    $x_min=$x[$ns];
	}
	if ($y_min>$y[$ns]) {
	    $y_min=$y[$ns];
	}
	if ($x_max<$x[$ns]) {
	    $x_max=$x[$ns];
	}
	if ($y_max<$y[$ns]) {
	    $y_max=$y[$ns];
	}
	if ($flux_min>$flux[$ns]) {
	    $flux_min=$flux[$ns];
	}
	if ($flux_max<$flux[$ns]) {
	    $flux_max=$flux[$ns];
	}
	
#	if ($ns>0) {
#	    if (($dpix>abs($x[$ns]-$x[$ns-1]))&&(abs($x[$ns]-$x[$ns-1])!=0)) {
#		$dpix=abs($x[$ns]-$x[$ns-1]);
#	    }
#	}
	$ns++;
#    }
}
close(FH);

$x_min=$x_min-$size;
$x_max=$x_max+$size;
$y_min=$y_min-$size;
$y_max=$y_max+$size;


$ns2=0;
open(FH,"<$slice2");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
	$id2[$ns2]=$data[0];
	$x2[$ns2]=$data[1];
	$y2[$ns2]=$data[2];
	$flux2[$ns2]=$data[3];
	$ns2++;
#    }
}
close(FH);


if ($#ARGV>=7) {
    $flux_min_plot=$ARGV[6];
    $flux_max_plot=$ARGV[7];
} else {
    $flux_min_plot=$flux_min;
    $flux_max_plot=$flux_max;
}
#print "$flux_min $flux_max\n";
for ($i=0;$i<$ns;$i++) {
    if ($flux_max_plot>$flux_min_plot) {
	$color[$i]=50+ceil((($flux[$i]-$flux_min_plot)/($flux_max_plot-$flux_min_plot))*50);
    } else {
	$color[$i]=1;
    }
#    $color[$i]=50+ceil(($flux[$i]/($flux_max_plot))*50);
    #$color[$i]=17+ceil(240*($flux[$i]/($flux_max_plot-$flux_min_plot)));
    #if ($flux[$i]>0) {
    #   print "$i $x[$i],$y[$i] $flux[$i] $color[$i]\n";
   #}
    $mask[$i]=0;
}

#$dev="/xs";
#$bright=1.0; $contrast=0.5;$table="smooth2"; $reverse=1; 
$bright=1.0; 
$contrast=0.5;
$table="idl5"; 
$reverse=1; 
if ($#ARGV==11) {
    $bright=$ARGV[8];
    $contrast=$ARGV[9];
    $table=$ARGV[10];
    $reverse=$ARGV[11];
}


pgbeg(0,$dev,1,1);
pgsubp(2,1);
pgscf(2.0);
pgsch(1.2);
($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
$nc=$pl->getdim(0);
for ($j=0;$j<$nc;$j++) {
    $l[$j] = $pl->slice($j)->sclr;
    $r[$j] = $pr->slice($j)->sclr;
    $g[$j] = $pg->slice($j)->sclr;
    $b[$j] = $pb->slice($j)->sclr;
}


$r[0]=1;  
$g[0]=1;  
$b[0]=1;  

pgscir(50,100);
pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

pgenv($x_max,$x_min,$y_min,$y_max,1,0);
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
for ($i=0;$i<$ns;$i++) {
    pgsci($color[$i]);

	if ($shape eq "C") {
	    pgsfs(1);
	    if ($flux[$i]>$flux_min_plot) {
		pgcirc($x[$i],$y[$i],$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
	    pgcirc($x[$i],$y[$i],$size/2);

	}
	if ($shape eq "R") {

	    pgsfs(1);
	    if ($flux[$i]>$flux_min_plot) {
		pgrect($x[$i]-$size/2,$x[$i]+$size/2,$y[$i]-$size/2,$y[$i]+$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
	    pgrect($x[$i]-$size/2,$x[$i]+$size/2,$y[$i]-$size/2,$y[$i]+$size/2);
	}
	    
}


pgsci(2);
if (($nmark>-1)&&($nmark<$ns)) {
    pgsch(4);
    pgslw(4);
    pgpoint(1,[$x[$nmark]],[$y[$nmark]],3);
}
pgsci(1);
pgwedg("RI", 0.1, 5, $flux_min_plot, $flux_max_plot, "");



# TWO

for ($i=0;$i<$ns2;$i++) {
    if ($flux_max_plot>$flux_min_plot) {
	$color2[$i]=50+ceil((($flux2[$i]-$flux_min_plot)/($flux_max_plot-$flux_min_plot))*50);
    } else {
	$color2[$i]=1;
    }
#    $color[$i]=50+ceil(($flux[$i]/($flux_max_plot))*50);
    #$color[$i]=17+ceil(240*($flux[$i]/($flux_max_plot-$flux_min_plot)));
    #if ($flux[$i]>0) {
    #   print "$i $x[$i],$y[$i] $flux[$i] $color[$i]\n";
   #}
    $mask[$i]=0;
}

pgenv($x_max,$x_min,$y_min,$y_max,1,0);
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
for ($i=0;$i<$ns2;$i++) {
    pgsci($color2[$i]);

	if ($shape eq "C") {
	    pgsfs(1);
	    if ($flux2[$i]>$flux_min_plot) {
		pgcirc($x2[$i],$y2[$i],$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
	    pgcirc($x2[$i],$y2[$i],$size/2);

	}
	if ($shape eq "R") {

	    pgsfs(1);
	    if ($flux[$i]>$flux_min_plot) {
		pgrect($x2[$i]-$size/2,$x2[$i]+$size/2,$y2[$i]-$size/2,$y2[$i]+$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
	    pgrect($x2[$i]-$size/2,$x2[$i]+$size/2,$y2[$i]-$size/2,$y2[$i]+$size/2);
	}
	    
}


pgsci(2);
if (($nmark>-1)&&($nmark<$ns)) {
    pgsch(4);
    pgslw(4);
    pgpoint(1,[$x2[$nmark]],[$y2[$nmark]],3);
}
pgsci(1);
pgwedg("RI", 0.1, 5, $flux_min_plot, $flux_max_plot, "");

pgclos();
pgend();




exit;
