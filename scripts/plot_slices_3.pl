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





if ($#ARGV<5) {
    print "USE: plot_slice.pl in_slice.txt in_slice2.txt in_slice3.txt DEVICE SHAPE SIZE [MARK] [MIN MAX] [BRIGHT CONTRAST LUT_TABLE INVERSE(0/1)]\n";
    exit;
}

$slice=$ARGV[0];
$slice2=$ARGV[1];
$slice3=$ARGV[2];
$dev=$ARGV[3];
$shape=$ARGV[4];
$size=$ARGV[5];

$nmark=-1;
if ($#ARGV>=6) {
    $nmark=$ARGV[6];
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
$x_c=0; $y_c=0; $sum_c=0;

open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	$flux[$ns]=$data[3];
    $l_flux[$ns]=log10($flux[$ns]);
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
	    $x_c=$x[$ns];
	    $y_c=$y[$ns];
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

$mean_flux=mean(@flux);
$r_min=1e12;
$r_max=-1e12;
$ns2=0;
open(FH,"<$slice2");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id2[$ns2]=$data[0];
    $x2[$ns2]=$data[1];
    $y2[$ns2]=$data[2];
    $flux2[$ns2]=$data[3];    
    $l_flux2[$ns2]=log10($flux2[$ns2]);
    $r[$ns2]=sqrt(($x[$ns2]-$x_c)**2+($y[$ns2]-$y_c)**2);
  #   if ($flux[$ns2]>$mean_flux) {
 #	$sum_c=$sum_c+$flux[$ns2];

  #  }



    $ns2++;
#    }
}
close(FH);


#if ($sum_c>0) {
#    $x_c=$x_c/$sum_c;
#    $y_c=$y_c/$sum_c;
#}


$ns3=0;
open(FH,"<$slice3");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
	$id3[$ns3]=$data[0];
	$x3[$ns3]=$data[1];
	$y3[$ns3]=$data[2];
	$flux3[$ns3]=$data[3];   
    $l_flux3[$ns3]=log10($flux3[$ns3]);

    

    if ($r_max<$r[$ns3]) {
	$r_max=$r[$ns3];
    }
    if ($r_min>$r[$ns3]) {
	$r_min=$r[$ns3];
    }


	$ns3++;
#    }
}
close(FH);

$mean_f1=mean(@flux);
$mean_f2=mean(@flux2);
$mean_f3=mean(@flux3);
$sigma_f3=sigma(@flux3);

$rat2=$mean_f2/$mean_f1;
$rat3=$mean_f3/$mean_f1;
$rat4=$sigma_f3/$mean_f1;
print "RAT = $rat2 $rat3 $rat4\n";

if ($#ARGV>=8) {
    $flux_min_plot=$ARGV[7];
    $flux_max_plot=$ARGV[8];
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
if ($#ARGV==12) {
    $bright=$ARGV[9];
    $contrast=$ARGV[10];
    $table=$ARGV[11];
    $reverse=$ARGV[12];
}


pgbeg(0,$dev,1,1);
pgsubp(2,2);
pgscf(2.0);
pgsch(1.3);
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
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","$slice");
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
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","$slice2");
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


# THREE

for ($i=0;$i<$ns3;$i++) {
    if ($flux_max_plot>$flux_min_plot) {
	$color3[$i]=50+ceil((($flux3[$i]-$flux_min_plot)/($flux_max_plot-$flux_min_plot))*50);
    } else {
	$color3[$i]=1;
    }
#    $color[$i]=50+ceil(($flux[$i]/($flux_max_plot))*50);
    #$color[$i]=17+ceil(340*($flux[$i]/($flux_max_plot-$flux_min_plot)));
    #if ($flux[$i]>0) {
    #   print "$i $x[$i],$y[$i] $flux[$i] $color[$i]\n";
   #}
    $mask[$i]=0;
}

pgenv($x_max,$x_min,$y_min,$y_max,1,0);
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","$slice3");
for ($i=0;$i<$ns3;$i++) {
    pgsci($color3[$i]);

	if ($shape eq "C") {
	    pgsfs(1);
	    if ($flux3[$i]>$flux_min_plot) {
		pgcirc($x3[$i],$y3[$i],$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
	    pgcirc($x3[$i],$y3[$i],$size/2);

	}
	if ($shape eq "R") {

	    pgsfs(1);
	    if ($flux[$i]>$flux_min_plot) {
		pgrect($x3[$i]-$size/2,$x3[$i]+$size/2,$y3[$i]-$size/2,$y3[$i]+$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
	    pgrect($x3[$i]-$size/2,$x3[$i]+$size/2,$y3[$i]-$size/2,$y3[$i]+$size/2);
	}
	    
}


pgsci(2);
if (($nmark>-1)&&($nmark<$ns)) {
    pgsch(4);
    pgslw(4);
    pgpoint(1,[$x3[$nmark]],[$y3[$nmark]],3);
}
pgsci(1);
pgwedg("RI", 0.1, 5, $flux_min_plot, $flux_max_plot, "");

$r_max=0.2*(abs($x_min)+abs($x_max));
$r_min=-0.1;
pgsci(1);
pgenv(-0.1*$flux_max,$flux_max,-0.1*$flux_max,$flux_max,1,0);
pgsch(5);
pgsci(1);
pgpoint($ns,\@flux,\@flux2,17);
pgsci(2);
pgpoint($ns,\@flux,\@flux3,22);

pgsci(1);
pgsch(1.3);
pglabel("Flux (data)","Flux (mod/res)","");

pgclos();
pgend();




exit;
