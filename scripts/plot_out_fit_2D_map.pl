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
    print "USE: plot_out_fit_2D_map.pl in_slice.txt mod_slice.txt DEVICE SHAPE SIZE [MARK] [MIN MAX] [BRIGHT CONTRAST LUT_TABLE INVERSE(0/1)]\n";
    exit;
}

$slice=$ARGV[0];
$mod_slice=$ARGV[1];
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

    if ($flux[$ns] eq "nan") {
	$flux[$ns]=0;
    }
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
	$ns++;
}
close(FH);

$ns=0;
open(FH,"<$mod_slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id[$ns]=$data[0];
    $x[$ns]=$data[1];
    $y[$ns]=$data[2];
    $flux_mod[$ns]=$data[3];
    if ($flux_mod[$ns] eq "nan") {
	$flux_mod[$ns]=0;
    }
    $flux_res[$ns]=$flux[$ns]-$flux_mos[$ns];
    $ns++;
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
	$color_mod[$i]=50+ceil((($flux_mod[$i]-$flux_min_plot)/($flux_max_plot-$flux_min_plot))*50);
	$color_res[$i]=50+ceil((($flux_res[$i]-$flux_min_plot)/($flux_max_plot-$flux_min_plot))*50);
    } else {
	$color[$i]=1;
    }
    $mask[$i]=0;
}

$bright=1.0; 
$contrast=0.5;
$table="real"; 
$reverse=1; 
if ($#ARGV==10) {
    $bright=$ARGV[7];
    $contrast=$ARGV[8];
    $table=$ARGV[9];
    $reverse=$ARGV[10];
}



pgbeg(0,$dev,1,1);
pgsubp(2,1);
pgscf(2.0);
pgsch(1.2);

$table="califaCT";
if ($table ne "califa") {
    if ($table eq "califa2") {
	while($cmap=<DATA>) {
	    chop($cmap);
	    @data=split(" ",$cmap);
	    $nc=$data[0];
	    $b[238-$nc]=$data[1];
	    $g[238-$nc]=$data[2];
	    $r[238-$nc]=$data[3];
	    $l[$nc]=$nc/238;
	}
	$bright=1;
	pgsitf(2);
    } else {
	if ($table eq "califaCT") {
	    $nc=254;
	    while($cmap=<DATA>) {
		chop($cmap);
		@data=split(" ",$cmap);
		$nc=$data[0];
#    $nc++;
		$r[$nc-1]=$data[1]/255;
		$g[$nc-1]=$data[2]/255;
		$b[$nc-1]=$data[3]/255;
		$l[$nc]=$nc/255;
#		print "";
#    $cmap=<DATA>
	    }
	    $bright=1; $contrast=0.5;
	    $r[0]=1.0;
	    $g[0]=1.0;
	    $b[0]=1.0;
	    $nc=254;

	} else {
	    ($pl,$pr,$pg,$pb)=lut_data($table,$reverse);
	    $nc=$pl->getdim(0);
	    for ($j=0;$j<$nc;$j++) {
		$l[$j] = $pl->slice($j)->sclr;
		$r[$j] = $pr->slice($j)->sclr;
		$g[$j] = $pg->slice($j)->sclr;
		$b[$j] = $pb->slice($j)->sclr;
	    }
	}
    }
} else {
    $nc=256;
    $R=-1.25;
    $S=0.5;
    $gamma=0.75;
    for ($j=0;$j<=256;$j++) {
	if ($reverse==0) {
	    $jj=$j;
	} else {
	    $jj=255-$j;
	}
	$l[$j] = $j*(1/256);
	$t=2*3.1416*($S/3+$R*$l[$j]);
	$a=(($l[$j])**$gamma)*(1-($l[$j])**$gamma)/2;
	$F=0.5+0.5*$l[$j];#exp(-((abs(($j-0)/256))**2));    
	$g[$jj] = $F*(($l[$j])**$gamma)-0.1486*$a*cos($t)+1.7828*$a*sin($t);
	$r[$jj] = $F*(($l[$j])**$gamma)-0.29223*$a*cos($t)-0.9065*$a*sin($t);
	$b[$jj] = $F*(($l[$j])**$gamma)+1.97294*$a*cos($t);
	}
}


for ($j=0;$j<$NW;$j++) {
    $r[$j]=1;
    $g[$j]=1;
    $b[$j]=1;
}

pgscir(50,100);
pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);
#pgclos();

#
# Original data
#
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
	    
#	    pgrect($x[$i]-$size/2,$x[$i]+$size/2,$y[$i]-$size/2,$y[$i]+$size/2);
	}
	    
}


pgsci(2);
if (($nmark>-1)&&($nmark<$ns)) {
#    pgsch(4);
#    pgslw(4);
    pgpoint(1,[$x[$nmark]],[$y[$nmark]],3);
}
pgsci(1);
pgwedg("RI", 0.1, 5, $flux_min_plot, $flux_max_plot, "");

#
# Model data
#
pgenv($x_max,$x_min,$y_min,$y_max,1,0);
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
for ($i=0;$i<$ns;$i++) {
    pgsci($color_mod[$i]);

	if ($shape eq "C") {
	    pgsfs(1);
	    if ($flux_mod[$i]>$flux_min_plot) {
		pgcirc($x[$i],$y[$i],$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
	    pgcirc($x[$i],$y[$i],$size/2);

	}
	if ($shape eq "R") {

	    pgsfs(1);
	    if ($flux_mod[$i]>$flux_min_plot) {
		pgrect($x[$i]-$size/2,$x[$i]+$size/2,$y[$i]-$size/2,$y[$i]+$size/2);
	    }
	    pgsfs(2);
	    pgsci(1);
	    
#	    pgrect($x[$i]-$size/2,$x[$i]+$size/2,$y[$i]-$size/2,$y[$i]+$size/2);
	}
	    
}


pgsci(2);
if (($nmark>-1)&&($nmark<$ns)) {
#    pgsch(4);
#    pgslw(4);
    pgpoint(1,[$x[$nmark]],[$y[$nmark]],3);
}
pgsci(1);
pgwedg("RI", 0.1, 5, $flux_min_plot, $flux_max_plot, "");

pgclos();
pgend();

print "****RESULTS****\n";
$call="cat out.fit_2D_map";
system($call);


exit;

__DATA__
      1.00000      146.543      0.00000      207.257
      2.00000      139.543      0.00000      206.257
      3.00000      135.314      0.00000      203.886
      4.00000      131.086      0.00000      201.514
      5.00000      126.857      0.00000      199.143
      6.00000      122.629      0.00000      196.771
      7.00000      118.400      0.00000      194.400
      8.00000      114.171      0.00000      192.029
      9.00000      109.943      0.00000      189.657
      10.0000      105.714      0.00000      187.286
      11.0000      101.486      0.00000      184.914
      12.0000      97.2571      0.00000      182.543
      13.0000      93.0286      0.00000      180.171
      14.0000      88.8000      0.00000      177.800
      15.0000      84.5714      0.00000      175.429
      16.0000      80.3429      0.00000      173.057
      17.0000      76.1143      0.00000      170.686
      18.0000      71.8857      0.00000      168.314
      19.0000      67.6572      0.00000      165.943
      20.0000      63.4286      0.00000      163.571
      21.0000      59.2000      0.00000      161.200
      22.0000      54.9714      0.00000      158.829
      23.0000      50.7429      0.00000      156.457
      24.0000      46.5143      0.00000      154.086
      25.0000      42.2857      0.00000      151.714
      26.0000      38.0571      0.00000      149.343
      27.0000      33.8286      0.00000      146.971
      28.0000      29.6000      0.00000      144.600
      29.0000      25.3714      0.00000      142.229
      30.0000      21.1429      0.00000      139.857
      31.0000      16.9143      0.00000      137.486
      32.0000      12.6857      0.00000      135.114
      33.0000      8.45715      0.00000      132.743
      34.0000      4.22858      0.00000      130.371
      35.0000      0.00000      0.00000      128.000
      36.0000      0.00000      3.47273      130.309
      37.0000      0.00000      6.94545      132.618
      38.0000      0.00000      10.4182      134.927
      39.0000      0.00000      13.8909      137.236
      40.0000      0.00000      17.3636      139.545
      41.0000      0.00000      20.8364      141.855
      42.0000      0.00000      24.3091      144.164
      43.0000      0.00000      27.7818      146.473
      44.0000      0.00000      31.2545      148.782
      45.0000      0.00000      34.7273      151.091
      46.0000      0.00000      38.2000      153.400
      47.0000      0.00000      41.6727      155.709
      48.0000      0.00000      45.1455      158.018
      49.0000      0.00000      48.6182      160.327
      50.0000      0.00000      52.0909      162.636
      51.0000      0.00000      55.5636      164.945
      52.0000      0.00000      59.0364      167.255
      53.0000      0.00000      62.5091      169.564
      54.0000      0.00000      65.9818      171.873
      55.0000      0.00000      69.4546      174.182
      56.0000      0.00000      72.9273      176.491
      57.0000      0.00000      76.4000      178.800
      58.0000      0.00000      79.8727      181.109
      59.0000      0.00000      83.3455      183.418
      60.0000      0.00000      86.8182      185.727
      61.0000      0.00000      90.2909      188.036
      62.0000      0.00000      93.7636      190.345
      63.0000      0.00000      97.2364      192.655
      64.0000      0.00000      100.709      194.964
      65.0000      0.00000      104.182      197.273
      66.0000      0.00000      107.655      199.582
      67.0000      0.00000      111.127      201.891
      68.0000      0.00000      114.600      204.200
      69.0000      0.00000      118.073      206.509
      70.0000      0.00000      121.545      208.818
      71.0000      0.00000      125.018      211.127
      72.0000      0.00000      128.491      213.436
      73.0000      0.00000      131.964      215.745
      74.0000      0.00000      135.436      218.055
      75.0000      0.00000      138.909      220.364
      76.0000      0.00000      142.382      222.673
      77.0000      0.00000      145.855      224.982
      78.0000      0.00000      149.327      227.291
      79.0000      0.00000      152.800      229.600
      80.0000      0.00000      156.273      231.909
      81.0000      0.00000      159.745      234.218
      82.0000      0.00000      163.218      236.527
      83.0000      0.00000      166.691      238.836
      84.0000      0.00000      170.164      241.145
      85.0000      0.00000      173.636      243.455
      86.0000      0.00000      177.109      245.764
      87.0000      0.00000      180.582      248.073
      88.0000      0.00000      184.055      250.382
      89.0000      0.00000      187.527      252.691
      90.0000      0.00000      191.000      255.000
      91.0000      1.57143      187.114      249.286
      92.0000      3.14286      183.229      243.571
      93.0000      4.71429      179.343      237.857
      94.0000      6.28571      175.457      232.143
      95.0000      7.85714      171.571      226.429
      96.0000      9.42857      167.686      220.714
      97.0000      11.0000      163.800      215.000
      98.0000      12.5714      159.914      209.286
      99.0000      14.1429      156.029      203.571
      100.000      15.7143      152.143      197.857
      101.000      17.2857      148.257      192.143
      102.000      18.8571      144.371      186.429
      103.000      20.4286      140.486      180.714
      104.000      22.0000      136.600      175.000
      105.000      23.5714      132.714      169.286
      106.000      25.1429      128.829      163.571
      107.000      26.7143      124.943      157.857
      108.000      28.2857      121.057      152.143
      109.000      29.8571      117.171      146.429
      110.000      31.4286      113.286      140.714
      111.000      33.0000      109.400      135.000
      112.000      34.5714      105.514      129.286
      113.000      36.1429      101.629      123.571
      114.000      37.7143      97.7429      117.857
      115.000      39.2857      93.8571      112.143
      116.000      40.8571      89.9714      106.429
      117.000      42.4286      86.0857      100.714
      118.000      44.0000      82.2000      95.0000
      119.000      45.5714      78.3143      89.2857
      120.000      47.1429      74.4286      83.5714
      121.000      48.7143      70.5428      77.8571
      122.000      50.2857      66.6571      72.1429
      123.000      51.8571      62.7714      66.4286
      124.000      53.4286      58.8857      60.7143
      125.000      55.0000      55.0000      55.0000
      126.000      59.7429      58.0000      59.7429
      127.000      64.4857      61.0000      64.4857
      128.000      69.2286      64.0000      69.2286
      129.000      73.9714      67.0000      73.9714
      130.000      78.7143      70.0000      78.7143
      131.000      83.4571      73.0000      83.4571
      132.000      88.2000      76.0000      88.2000
      133.000      92.9429      79.0000      92.9429
      134.000      97.6857      82.0000      97.6857
      135.000      102.429      85.0000      102.429
      136.000      107.171      88.0000      107.171
      137.000      111.914      91.0000      111.914
      138.000      116.657      94.0000      116.657
      139.000      121.400      97.0000      121.400
      140.000      126.143      100.000      126.143
      141.000      130.886      103.000      130.886
      142.000      135.629      106.000      135.629
      143.000      140.371      109.000      140.371
      144.000      145.114      112.000      145.114
      145.000      149.857      115.000      149.857
      146.000      154.600      118.000      154.600
      147.000      159.343      121.000      159.343
      148.000      164.086      124.000      164.086
      149.000      168.829      127.000      168.829
      150.000      173.571      130.000      173.571
      151.000      178.314      133.000      178.314
      152.000      183.057      136.000      183.057
      153.000      187.800      139.000      187.800
      154.000      192.543      142.000      192.543
      155.000      197.286      145.000      197.286
      156.000      202.029      148.000      202.029
      157.000      206.771      151.000      206.771
      158.000      211.514      154.000      211.514
      159.000      216.257      157.000      216.257
      160.000      221.000      160.000      221.000
      161.000      221.567      157.333      217.317
      162.000      222.133      154.667      213.633
      163.000      222.700      152.000      209.950
      164.000      223.267      149.333      206.267
      165.000      223.833      146.667      202.583
      166.000      224.400      144.000      198.900
      167.000      224.967      141.333      195.217
      168.000      225.533      138.667      191.533
      169.000      226.100      136.000      187.850
      170.000      226.667      133.333      184.167
      171.000      227.233      130.667      180.483
      172.000      227.800      128.000      176.800
      173.000      228.367      125.333      173.117
      174.000      228.933      122.667      169.433
      175.000      229.500      120.000      165.750
      176.000      230.067      117.333      162.067
      177.000      230.633      114.667      158.383
      178.000      231.200      112.000      154.700
      179.000      231.767      109.333      151.017
      180.000      232.333      106.667      147.333
      181.000      232.900      104.000      143.650
      182.000      233.467      101.333      139.967
      183.000      234.033      98.6667      136.283
      184.000      234.600      96.0000      132.600
      185.000      235.167      93.3333      128.917
      186.000      235.733      90.6667      125.233
      187.000      236.300      88.0000      121.550
      188.000      236.867      85.3333      117.867
      189.000      237.433      82.6667      114.183
      190.000      238.000      80.0000      110.500
      191.000      238.567      77.3333      106.817
      192.000      239.133      74.6667      103.133
      193.000      239.700      72.0000      99.4500
      194.000      240.267      69.3333      95.7667
      195.000      240.833      66.6667      92.0833
      196.000      241.400      64.0000      88.4000
      197.000      241.967      61.3333      84.7167
      198.000      242.533      58.6667      81.0333
      199.000      243.100      56.0000      77.3500
      200.000      243.667      53.3333      73.6667
      201.000      244.233      50.6667      69.9833
      202.000      244.800      48.0000      66.3000
      203.000      245.367      45.3333      62.6167
      204.000      245.933      42.6667      58.9333
      205.000      246.500      40.0000      55.2500
      206.000      247.067      37.3333      51.5667
      207.000      247.633      34.6667      47.8833
      208.000      248.200      32.0000      44.2000
      209.000      248.767      29.3333      40.5167
      210.000      249.333      26.6667      36.8333
      211.000      249.900      24.0000      33.1500
      212.000      250.467      21.3333      29.4667
      213.000      251.033      18.6667      25.7833
      214.000      251.600      16.0000      22.1000
      215.000      252.167      13.3333      18.4167
      216.000      252.733      10.6667      14.7333
      217.000      253.300      8.00000      11.0500
      218.000      253.867      5.33333      7.36666
      219.000      254.433      2.66667      3.68332
      220.000      255.000      0.00000      0.00000
      221.000      255.000      4.71429      0.00000
      222.000      255.000      9.42857      0.00000
      223.000      255.000      14.1429      0.00000
      224.000      255.000      18.8571      0.00000
      225.000      255.000      23.5714      0.00000
      226.000      255.000      28.2857      0.00000
      227.000      255.000      33.0000      0.00000
      228.000      255.000      37.7143      0.00000
      229.000      255.000      42.4286      0.00000
      230.000      255.000      47.1429      0.00000
      231.000      255.000      51.8571      0.00000
      232.000      255.000      56.5714      0.00000
      233.000      255.000      61.2857      0.00000
      234.000      255.000      66.0000      0.00000
      235.000      255.000      70.7143      0.00000
      236.000      255.000      75.4286      0.00000
      237.000      255.000      80.1429      0.00000
      238.000      255.000      84.8571      0.00000
      239.000      255.000      89.5714      0.00000
      240.000      255.000      94.2857      0.00000
      241.000      255.000      99.0000      0.00000
      242.000      255.000      103.714      0.00000
      243.000      255.000      108.429      0.00000
      244.000      255.000      113.143      0.00000
      245.000      255.000      117.857      0.00000
      246.000      255.000      122.571      0.00000
      247.000      255.000      127.286      0.00000
      248.000      255.000      132.000      0.00000
      249.000      255.000      136.714      0.00000
      250.000      255.000      141.429      0.00000
      251.000      255.000      146.143      0.00000
      252.000      255.000      150.857      0.00000
      253.000      255.000      155.571      0.00000
      254.000      255.000      160.286      0.00000
      255.000      255.000      165.000      0.00000
      256.000      255.000      165.000      0.00000

