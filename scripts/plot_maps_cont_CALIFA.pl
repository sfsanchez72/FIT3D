#!/usr/bin/env perl
use Statistics::OLS;
#use Math::FFT;
#use Math::Stat;
#use Math::Spline qw(spline linsearch binsearch);
#use Math::Derivative qw(Derivative2);
#use Math::Approx;
use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Core;
use PDL::Graphics::LUT;
use Carp;



$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "my.pl";

if ($#ARGV<13) {
    print "USE: maps_orion.pl map.fits min max bright contrast Label factor dev SYSTEMIC_VEL MAP_CONT NLEVELS MIN STEP NWHITE[49] [ROT]\n";
    exit;
}
$mapfile=$ARGV[0];
$min=$ARGV[1];
$max=$ARGV[2];
$bright=$ARGV[3];
$contrast=$ARGV[4];
$label=$ARGV[5];
$factor=$ARGV[6];
$dev=$ARGV[7];
$vw=$ARGV[8];
$map_cont=$ARGV[9];
$nlevels=$ARGV[10];
$min_cont=$ARGV[11];
$steps=$ARGV[12];
$NW=$ARGV[13];

if ($dev !~ "TPNG") { 
    $ENV{PGPLOT_FOREGROUND} = "black";
    $ENV{PGPLOT_BACKGROUND} = "white";
} else {
    $ENV{PGPLOT_FOREGROUND} = "white";
    $ENV{PGPLOT_BACKGROUND} = "white";
}



$rot=0;
if ($#ARGV==14) {
     $rot=$ARGV[14];
}

$color_cont=0;

#$pdl_mask=rfits("M74_mos.mask.fits");
#$pdl_mask=rfits("mask_sc.fits");
#$pdl_map=rfits($mapfile);
#$pdl_cont=rfits($map_cont);



$pdl_map_tmp=rfits($mapfile);
$pdl_cont_tmp=rfits($map_cont);

$is_zero=sum(abs($pdl_map_tmp));



if ($rot==0) {
    $pdl_map=$pdl_map_tmp;
    $pdl_cont=$pdl_cont_tmp;
} else {
#    $pdl_map=$pdl_map_tmp->rot2d($rot,0,1);
#    $pdl_cont=$pdl_cont_tmp->rot2d($rot,0,1);
    $pdl_map=transpose($pdl_map_tmp);
    $pdl_cont=transpose($pdl_cont_tmp);
}



#$crval1=$pdl_map->hdr->{CRVAL1};
#$cdelt1=$pdl_map->hdr->{CDELT1};
#$crval2=$pdl_map->hdr->{CRVAL2};
#$cdelt2=$pdl_map->hdr->{CDELT2};
#print "$crval1 $cdelt1 $crval2 $cdelt2\n";
($nx,$ny)=$pdl_map->dims;
$h=$pdl_map->gethdr;
$crval1=-138.839996337891;
$cdelt1=1;
$crval2=-210.726806640625;
$cdelt2=1;


$pdl_map=$pdl_map*$factor;


$NNx=$nx-1;
$NNy=$ny-1;
$rpdl_map=$pdl_map->slice("$NNx:0,$NNy:0");
@map=list($rpdl_map);
$rpdl_cont=$pdl_cont->slice("$NNx:0,$NNy:0");
@cont=list($rpdl_cont);
$crval1=-30;
$crval2=-40;
$x_max=$crval1;
$x_min=$crval1+$cdelt1*$nx;
$y_min=$crval2;
$y_max=$crval2+$cdelt2*$ny;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	if ($map[$i+$j*$nx]>0) {
	    $map[$i+$j*$nx]=$map[$i+$j*$nx]-$vw;
	}

	if ($map[$i+$j*$nx]>1e12) {
	    $map[$i+$j*$nx]=0;
	}
	if ($map[$i+$j*$nx] eq "nan") {
	    $map[$i+$j*$nx]=0;
	}


	if ($map[$i+$j*$nx]!=(1*$map[$i+$j*$nx])) {
	    $map[$i+$j*$nx]=0;
	}
	if ($map[$i+$j*$nx]==$vw) {
	    $map[$i+$j*$nx]=-1e12;
	}
#	print "$map[$i+$j*$nx]\n";

    }
}


@tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
#@tr=($crval1,$cdelt1,0,$crval2,0,$cdelt2);
#@tr=($x_max,-$cdelt1,0,$crval2,0,$cdelt2);



#@tr=(0,1,0,0,0,1);
#$x_min=0;
#$x_max=$nx;
#$y_min=0;
#$y_max=$ny;

pgbeg(0,$dev,1,1);
pgscf(2.0);
pgsch(1.2);


$nc0=255;
while($cmap=<DATA>) {
    chop($cmap);
    @data=split(" ",$cmap);
    $nc=int($data[0]);
    $k=256-$nc-1;
    $r[$k]=$data[1]/255;
    $g[$k]=$data[2]/255;
    $b[$k]=$data[3]/255;

    $k=$nc-1;
    $rr[$k]=$data[1]/255;
    $gg[$k]=$data[2]/255;
    $bb[$k]=$data[3]/255;

    $l[$nc-1]=$nc/255;
#    print "$nc $k $r[$k] $g[$k] $b[$k]\n";
#    $cmap=<DATA>;
}
#$bright=1; $contrast=0.5;
$rr[0]=1.0;
$gg[0]=1.0;
$bb[0]=1.0;
$r[0]=1.0;
$g[0]=1.0;
$b[0]=1.0;

pgscir(50,128+50); $nc=254;
#pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);

pgctab(\@l,\@rr,\@gg,\@bb,$nc,$bright,$contrast);

#pgenv($x_min,$x_max,$y_min,$y_max,1,0);
pgslw(2);
pgenv($x_min,$x_max,$y_min,$y_max,1,0);
#pgpap(6.0,1.0);

#pgsch(0.7);           # Set character height
#pgsch(1.2);

if ($dev !~ "TPNG") { 
pglabel("\\gD RA (arcsec)","\\gD DEC (arcsec)","");
}
if ($is_zero!=0) {
    pgimag(\@map,$nx,$ny,1,$nx,1,$ny,$min,$max,\@tr);
}
pgsci(0);
pgsch(1);
#pgsch(1.1);
@tmp;
$nt=0;
for ($j=0;$j<$ny;$j++) {
    for ($i=0;$i<$nx;$i++) {
	$val=$map[$i+$j*$nx];
#	if ($map[$i+$j*$nx]<=0) {
	if ($val<=0) {
	    $X=$crval1+$cdelt1*$i+1;
	    $Y=$crval2+$cdelt2*$j+1;
	#    print "$i,$j $X,$Y MAP=$map[$i+$j*$nx]\n";
#	    pgpoint(1,[$X],[$Y],16);
	} else {
	    if ($val>6) {
		$tmp[$nt]=$val;
#		print "$tmp[$nt] $nt\n";
		$nt++;
	    }
	}
    }
}
$mean=mean(@tmp);
$sigma=sigma(@tmp);
$sum=$nt*$mean;
print "FLUX=$sum\n";
print "MEAN=$mean+-$sigma\n";


pgsch(3);
pgsci(3);
pgslw(3);
#pgpoint1,[3.65],[-3.13],3);
#pgpoint(1,0,0,3);
pgsch(1.2);
pgslw(1);


pgsci(1);
pgslw(3);
#$nlevels=10;
#$min_cont=$ARGV[11];
#$steps=$ARGV[12];
for ($i=0;$i<$nlevels;$i++) {    
    $levels[$i]=$min_cont+$steps*($i**1.5);
#    print "$i $levels[$i]\n";
}

pgsci(1);
pgslw(4);

if ($dev !~ "TPNG") { 
    pgsci(1);
    pgslw(4);
} else { 
    pgsci(14);
    pgslw(8);
}

if ($is_zero!=0) {
  pgsci(1);
  pgslw(4);
    pgcont(\@cont,$nx,$ny,1,$nx,1,$ny,\@levels,$nlevels,\@tr);
  
}
pgslw(1);
pgsci(1);

if ($dev !~ "TPNG") { 
    pgsci(1);
    pgsch(1.4);
    pgptxt(-15,28,0,0.5,$label);
    pgsci(1);


pgsch(1.1);
#pgptxt(-85,160,0,0,$label);
#pgptxt(-110,160,0,0,$label);

for ($j=0;$j<256;$j++) {
#    if (abs($j-128)>2) {
	$F=0.85/exp(-(($j-128)/512)**4);    
#    } else {
#	$F=1;       
#    }
    $rr[$j]=$r[$j]*$F;
    $gg[$j]=$g[$j]*$F;
    $bb[$j]=$b[$j]*$F;
    
}
#pgctab(\@l,\@rr,\@gg,\@bb,$nc,$c,0.5);
pgslw(2);
pgwedg("RI",0.0,4.5,$min,$max,"");
pgbox("SBC",0,0,"SBC",0,0);
}




pgclos();
pgend();



exit;


__DATA__
      1.00000      0.00000      5.64000      130.080
      2.00000      0.00000      7.64000      133.080
      3.00000      0.00000      11.4600      135.620
      4.00000      0.00000      15.2800      138.160
      5.00000      0.00000      19.1000      140.700
      6.00000      0.00000      22.9200      143.240
      7.00000      0.00000      26.7400      145.780
      8.00000      0.00000      30.5600      148.320
      9.00000      0.00000      34.3800      150.860
      10.0000      0.00000      38.2000      153.400
      11.0000      0.00000      42.0200      155.940
      12.0000      0.00000      45.8400      158.480
      13.0000      0.00000      49.6600      161.020
      14.0000      0.00000      53.4800      163.560
      15.0000      0.00000      57.3000      166.100
      16.0000      0.00000      61.1200      168.640
      17.0000      0.00000      64.9400      171.180
      18.0000      0.00000      68.7600      173.720
      19.0000      0.00000      72.5800      176.260
      20.0000      0.00000      76.4000      178.800
      21.0000      0.00000      80.2200      181.340
      22.0000      0.00000      84.0400      183.880
      23.0000      0.00000      87.8600      186.420
      24.0000      0.00000      91.6800      188.960
      25.0000      0.00000      95.5000      191.500
      26.0000      0.00000      99.3200      194.040
      27.0000      0.00000      103.140      196.580
      28.0000      0.00000      106.960      199.120
      29.0000      0.00000      110.780      201.660
      30.0000      0.00000      114.600      204.200
      31.0000      0.00000      118.420      206.740
      32.0000      0.00000      122.240      209.280
      33.0000      0.00000      126.060      211.820
      34.0000      0.00000      129.880      214.360
      35.0000      0.00000      133.700      216.900
      36.0000      0.00000      137.520      219.440
      37.0000      0.00000      141.340      221.980
      38.0000      0.00000      145.160      224.520
      39.0000      0.00000      148.980      227.060
      40.0000      0.00000      152.800      229.600
      41.0000      0.00000      156.620      232.140
      42.0000      0.00000      160.440      234.680
      43.0000      0.00000      164.260      237.220
      44.0000      0.00000      168.080      239.760
      45.0000      0.00000      171.900      242.300
      46.0000      0.00000      175.720      244.840
      47.0000      0.00000      179.540      247.380
      48.0000      0.00000      183.360      249.920
      49.0000      0.00000      187.180      252.460
      50.0000      0.00000      191.000      255.000
      51.0000      5.10000      187.180      249.900
      52.0000      10.2000      183.360      244.800
      53.0000      15.3000      179.540      239.700
      54.0000      20.4000      175.720      234.600
      55.0000      25.5000      171.900      229.500
      56.0000      30.6000      168.080      224.400
      57.0000      35.7000      164.260      219.300
      58.0000      40.8000      160.440      214.200
      59.0000      45.9000      156.620      209.100
      60.0000      51.0000      152.800      204.000
      61.0000      56.1000      148.980      198.900
      62.0000      61.2000      145.160      193.800
      63.0000      66.3000      141.340      188.700
      64.0000      71.4000      137.520      183.600
      65.0000      76.5000      133.700      178.500
      66.0000      81.6000      129.880      173.400
      67.0000      86.7000      126.060      168.300
      68.0000      91.8000      122.240      163.200
      69.0000      96.9000      118.420      158.100
      70.0000      102.000      114.600      153.000
      71.0000      107.100      110.780      147.900
      72.0000      112.200      106.960      142.800
      73.0000      117.300      103.140      137.700
      74.0000      122.400      99.3200      132.600
      75.0000      127.500      95.5000      127.500
      76.0000      132.600      91.6800      122.400
      77.0000      137.700      87.8600      117.300
      78.0000      142.800      84.0400      112.200
      79.0000      147.900      80.2200      107.100
      80.0000      153.000      76.4000      102.000
      81.0000      158.100      72.5800      96.9000
      82.0000      163.200      68.7600      91.8000
      83.0000      168.300      64.9400      86.7000
      84.0000      173.400      61.1200      81.6000
      85.0000      178.500      57.3000      76.5000
      86.0000      183.600      53.4800      71.4000
      87.0000      188.700      49.6600      66.3000
      88.0000      193.800      45.8400      61.2000
      89.0000      198.900      42.0200      56.1000
      90.0000      204.000      38.2000      51.0000
      91.0000      209.100      34.3800      45.9000
      92.0000      214.200      30.5600      40.8000
      93.0000      219.300      26.7400      35.7000
      94.0000      224.400      22.9200      30.6000
      95.0000      229.500      19.1000      25.5000
      96.0000      234.600      15.2800      20.4000
      97.0000      239.700      11.4600      15.3000
      98.0000      244.800      7.64001      10.2000
      99.0000      249.900      3.82000      5.10000
      100.000      255.000      0.00000      0.00000
      101.000      255.000      3.30000      0.00000
      102.000      255.000      6.60000      0.00000
      103.000      255.000      9.90000      0.00000
      104.000      255.000      13.2000      0.00000
      105.000      255.000      16.5000      0.00000
      106.000      255.000      19.8000      0.00000
      107.000      255.000      23.1000      0.00000
      108.000      255.000      26.4000      0.00000
      109.000      255.000      29.7000      0.00000
      110.000      255.000      33.0000      0.00000
      111.000      255.000      36.3000      0.00000
      112.000      255.000      39.6000      0.00000
      113.000      255.000      42.9000      0.00000
      114.000      255.000      46.2000      0.00000
      115.000      255.000      49.5000      0.00000
      116.000      255.000      52.8000      0.00000
      117.000      255.000      56.1000      0.00000
      118.000      255.000      59.4000      0.00000
      119.000      255.000      62.7000      0.00000
      120.000      255.000      66.0000      0.00000
      121.000      255.000      69.3000      0.00000
      122.000      255.000      72.6000      0.00000
      123.000      255.000      75.9000      0.00000
      124.000      255.000      79.2000      0.00000
      125.000      255.000      82.5000      0.00000
      126.000      255.000      85.8000      0.00000
      127.000      255.000      89.1000      0.00000
      128.000      255.000      92.4000      0.00000
      129.000      255.000      95.7000      0.00000
      130.000      255.000      99.0000      0.00000
      131.000      255.000      102.300      0.00000
      132.000      255.000      105.600      0.00000
      133.000      255.000      108.900      0.00000
      134.000      255.000      112.200      0.00000
      135.000      255.000      115.500      0.00000
      136.000      255.000      118.800      0.00000
      137.000      255.000      122.100      0.00000
      138.000      255.000      125.400      0.00000
      139.000      255.000      128.700      0.00000
      140.000      255.000      132.000      0.00000
      141.000      255.000      135.300      0.00000
      142.000      255.000      138.600      0.00000
      143.000      255.000      141.900      0.00000
      144.000      255.000      145.200      0.00000
      145.000      255.000      148.500      0.00000
      146.000      255.000      151.800      0.00000
      147.000      255.000      155.100      0.00000
      148.000      255.000      158.400      0.00000
      149.000      255.000      161.700      0.00000
      150.000      255.000      165.000      0.00000
      151.000      251.000      162.800      1.10000
      152.000      247.000      160.600      2.20000
      153.000      243.000      158.400      3.30000
      154.000      239.000      156.200      4.40000
      155.000      235.000      154.000      5.50000
      156.000      231.000      151.800      6.60000
      157.000      227.000      149.600      7.70000
      158.000      223.000      147.400      8.80000
      159.000      219.000      145.200      9.90000
      160.000      215.000      143.000      11.0000
      161.000      211.000      140.800      12.1000
      162.000      207.000      138.600      13.2000
      163.000      203.000      136.400      14.3000
      164.000      199.000      134.200      15.4000
      165.000      195.000      132.000      16.5000
      166.000      191.000      129.800      17.6000
      167.000      187.000      127.600      18.7000
      168.000      183.000      125.400      19.8000
      169.000      179.000      123.200      20.9000
      170.000      175.000      121.000      22.0000
      171.000      171.000      118.800      23.1000
      172.000      167.000      116.600      24.2000
      173.000      163.000      114.400      25.3000
      174.000      159.000      112.200      26.4000
      175.000      155.000      110.000      27.5000
      176.000      151.000      107.800      28.6000
      177.000      147.000      105.600      29.7000
      178.000      143.000      103.400      30.8000
      179.000      139.000      101.200      31.9000
      180.000      135.000      99.0000      33.0000
      181.000      131.000      96.8000      34.1000
      182.000      127.000      94.6000      35.2000
      183.000      123.000      92.4000      36.3000
      184.000      119.000      90.2000      37.4000
      185.000      115.000      88.0000      38.5000
      186.000      111.000      85.8000      39.6000
      187.000      107.000      83.6000      40.7000
      188.000      103.000      81.4000      41.8000
      189.000      99.0000      79.2000      42.9000
      190.000      95.0000      77.0000      44.0000
      191.000      91.0000      74.8000      45.1000
      192.000      87.0000      72.6000      46.2000
      193.000      83.0000      70.4000      47.3000
      194.000      79.0000      68.2000      48.4000
      195.000      75.0000      66.0000      49.5000
      196.000      71.0000      63.8000      50.6000
      197.000      67.0000      61.6000      51.7000
      198.000      63.0000      59.4000      52.8000
      199.000      59.0000      57.2000      53.9000
      200.000      55.0000      55.0000      55.0000
      201.000      58.0182      56.9091      58.0182
      202.000      61.0364      58.8182      61.0364
      203.000      64.0545      60.7273      64.0545
      204.000      67.0727      62.6364      67.0727
      205.000      70.0909      64.5455      70.0909
      206.000      73.1091      66.4545      73.1091
      207.000      76.1273      68.3636      76.1273
      208.000      79.1455      70.2727      79.1455
      209.000      82.1636      72.1818      82.1636
      210.000      85.1818      74.0909      85.1818
      211.000      88.2000      76.0000      88.2000
      212.000      91.2182      77.9091      91.2182
      213.000      94.2364      79.8182      94.2364
      214.000      97.2545      81.7273      97.2545
      215.000      100.273      83.6364      100.273
      216.000      103.291      85.5455      103.291
      217.000      106.309      87.4546      106.309
      218.000      109.327      89.3636      109.327
      219.000      112.345      91.2727      112.345
      220.000      115.364      93.1818      115.364
      221.000      118.382      95.0909      118.382
      222.000      121.400      97.0000      121.400
      223.000      124.418      98.9091      124.418
      224.000      127.436      100.818      127.436
      225.000      130.455      102.727      130.455
      226.000      133.473      104.636      133.473
      227.000      136.491      106.545      136.491
      228.000      139.509      108.455      139.509
      229.000      142.527      110.364      142.527
      230.000      145.545      112.273      145.545
      231.000      148.564      114.182      148.564
      232.000      151.582      116.091      151.582
      233.000      154.600      118.000      154.600
      234.000      157.618      119.909      157.618
      235.000      160.636      121.818      160.636
      236.000      163.655      123.727      163.655
      237.000      166.673      125.636      166.673
      238.000      169.691      127.545      169.691
      239.000      172.709      129.455      172.709
      240.000      175.727      131.364      175.727
      241.000      178.745      133.273      178.745
      242.000      181.764      135.182      181.764
      243.000      184.782      137.091      184.782
      244.000      187.800      139.000      187.800
      245.000      190.818      140.909      190.818
      246.000      193.836      142.818      193.836
      247.000      196.855      144.727      196.855
      248.000      199.873      146.636      199.873
      249.000      202.891      148.545      202.891
      250.000      205.909      150.455      205.909
      251.000      208.927      152.364      208.927
      252.000      211.945      154.273      211.945
      253.000      214.964      156.182      214.964
      254.000      217.982      158.091      217.982
      255.000      221.000      160.000      221.000
      256.000      255.000      165.000      0.00000
