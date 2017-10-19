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


if ($#ARGV<5) {
    print "USE: get_res_plot.pl ASCII_TABLE.txt CRVAL CDELT LABEL1 LABEL2 DEV [MIN MAX] [YMIN YMAX] [NMAX] [NPOLY] [FINAL_RESOLUTION] [OUTPUTFILE]\n";
    exit;
}

$infile=$ARGV[0];
$nc1=$ARGV[1];
$nc2=$ARGV[2];
$label1=$ARGV[3];
$label2=$ARGV[4];
$dev=$ARGV[5];

$crval=$nc1;
$cdelt=$nc2;
$npoly=7;
$nmax=1970;
$final_res=3.3;
$output="get_res_plot.fits";


$N=0;
open(NL,"wc -l $infile |");
$NL=<NL>;
chop($NL);
($N,$name)=split(" ",$NL);
close(NL);


$n=0;
$K=0;
open(FH,"<$infile");
for ($i=0;$i<$nc1;$i++) {
    $line=<FH>;
    chop($line);
    $nc2=$line;
    $nc1=$N/($line+1);
    for ($j=0;$j<$nc2;$j++) {
	$line=<FH>;
	chop($line);
	if ($line !~ "#") {
	    @data=split(" ",$line);
	    $x[$n]=$data[1];	
	    $y[$n]=$data[5];
	    $delta[$n]=$y[$n];
	    $n++;
	    
	}
    }
    $slice=pdl(@x);
    ($mean,$rms,$median,$min,$max) = stats($slice);
    $x_min=$min-0.2*($median-$min);
    $x_max=$max+0.2*($max-$median);
    $label12=$label1." ".$i."/".$nc1;
    $slice=pdl(@y);
    ($mean,$rms,$median,$min,$max) = stats($slice);
    $y_min=0.99*$min;
    $y_max=1.01*$max;
    $NN=$n/int($nbin/6+1);
    
    
    $slice=pdl(@delta);
    ($mean,$rms,$median,$min,$max) = stats($slice);
    
#print "Y-X = $mean $rms $median $min $max\n";

    if ($#ARGV==7) {
	$x_min=$ARGV[6];
	$x_max=$ARGV[7];
    }
    
    if ($#ARGV==9) {
	$x_min=$ARGV[6];
	$x_max=$ARGV[7];
	$y_min=$ARGV[8];
	$y_max=$ARGV[9];
    }
    

$fix=0;

if ($#ARGV==10) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $nmax=$ARGV[10];
}

if ($#ARGV==11) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $nmax=$ARGV[10];
    $npoly=$ARGV[11];
}

if ($#ARGV==12) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $nmax=$ARGV[10];
    $npoly=$ARGV[11];
    $final_res=$ARGV[12];
}

if ($#ARGV==13) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $nmax=$ARGV[10];
    $npoly=$ARGV[11];
    $final_res=$ARGV[12];
    $output=$ARGV[13];
}



    if ($K==0) {
	$pdl_out=zeroes($nmax,$nc1);
    }

    for ($ii=0;$ii<$nmax;$ii++) {
	$wave[$ii]=$crval+$cdelt*$ii;
	$flux[$ii]=0;
    }

    pgbegin(0,$dev,1,1);
    pgsfs(1.2);
    pgscf(2);             # Set character font
    pgslw(2);             # Set line width
    pgsch(1.4);           # Set character height
    pgenv($x_min,$x_max,$y_min,$y_max,$fix,0);
    pglabel("$label12","$label2","");
    pgsch(2.5);
    $k=0;
#print "$n\n @x\n @y\n";
    for ($ii=0;$ii<$n;$ii++) {
	$X=$x[$ii];
	$Y=$y[$ii];
	$D[$ii]=$Y;
	if ($Y!=0) {
	    $R[$k]=$X/$Y;
	    $k++;
	}
#    print "$i $X $Y\n";
	pgsci(4);
	pgpoint(1,[$X],[$Y],17);
	pgsci(1);
	pgpoint(1,[$X],[$Y],22);
    }
#pgpoint($n,\@x,\@y,16);
    

    $med=median(@D);
    $sig=sigma(@D);
    ($min,$max)=minmax(@D);
#    print "Y =$med +- $sig $min $max\n";
    
    $nk=0;
    $xp[$nk]=$wave[0];
    $yp[$nk]=$med;
    $nk++;
    for ($ii=0;$ii<$n;$ii++) {
	if (abs($y[$ii]-$med)<2*$sig) {
	    $yp[$nk]=$y[$ii];
	    $xp[$nk]=$x[$ii];	    
	    if (($yp[$nk]<1.25)&&($yp[$nk]>0.7)&&($xp[$nk]>3650)) {
		$nk++;
	    }
	}
    }
    pgsci(8);
    pgpoint($nk,\@xp,\@yp,2);
    pgsci(1);
#    $xp[$nk]=$wave[$nmax-1];
#    $yp[$nk]=$med;
    ($s_x,$coeff) = fitpoly1d(pdl(@xp),pdl(@yp),$npoly);
    @c=list($coeff);
    for ($ii=0;$ii<$nmax;$ii++) {
	$out_spec[$ii]=0;
	for ($k=0;$k<$npoly;$k++) {
            $out_spec[$ii]=$out_spec[$ii]+$c[$k]*(($wave[$ii])**$k);	   
        }
	$val=$out_spec[$ii];
	if (($val**2)>($final_res**2)) {
	    $final_val=0;
	} else {
	    $final_val=sqrt($final_res**2-$val**2);
	}
#	print "$ii $wave[$ii] $final_val\n";
	set($pdl_out,$ii,$i,$final_val);
    }
    
    pgsci(2);
    pgline($nmax,\@wave,\@out_spec);
    pgsci(1);
#    $label22=$label2." ".$med."+-".$sig;
    pgsch(1.2);
    pgptxt(1.05*$x_min,0.9*$y_max,0,0,"$med+-$sig");
#    pgenv($x_min,$x_max,$y_min,$y_max,$fix,0);





#    print "Press Enter"; <stdin>;
    $n=0;
    $K++;
}
    pgclose;
    pgend;
    close(FH);

$h = {NAXIS=>2, NAXIS1=>$nmax, NAXIS=>$nc1, CRPIX => 1, CRVAL1 => $crval, CDELT1 => $cdelt, COMMENT=>"Sample FITS-style header"};

$pdl_out->sethdr($h);

$pdl_out->wfits($output);

    exit;





