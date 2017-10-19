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
    print "USE: table_plot.pl ASCII_TABLE.txt N.COLUMN1(0...N)  N.COLUMN2(0...N) LABEL1 LABEL2 DEV [MIN MAX] [YMIN YMAX] [FIX] [ONE-to-ONE_LINE 0/1]\n";
    exit;
}

$infile=$ARGV[0];
$nc1=$ARGV[1];
$nc2=$ARGV[2];
$label1=$ARGV[3];
$label2=$ARGV[4];
$dev=$ARGV[5];
$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$x[$n]=$data[$nc1];	
	$y[$n]=$data[$nc2];
	$delta[$n]=$y[$n]-$x[$n];
	if ($line !~ "LINE") {
	    $n++;
	}
    }
}
close(FH);
$slice=pdl(@x);
($mean,$rms,$median,$min,$max) = stats($slice);
$median=median(@x);
$x_min=$min-0.2*($median-$min);
$x_max=$max+0.2*($max-$median);
if ($median>2000000) {
    $label1=$label1."-".$median;
    for ($i=0;$i<$n;$i++) {
	$x[$i]=$x[$i]-$median;
    }
    $x_min=$x_min-$median;
    $x_max=$x_max-$median;
}
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
    $fix=$ARGV[10];
}

$one=0;
if ($#ARGV==11) {
    $x_min=$ARGV[6];
    $x_max=$ARGV[7];
    $y_min=$ARGV[8];
    $y_max=$ARGV[9];
    $fix=$ARGV[10];
    $one=$ARGV[11];
}


pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv($x_min,$x_max,$y_min,$y_max,$fix,0);
pglabel("$label1","$label2","");
if ($one==1) {
    pgsci(2);
    pgsls(2);
    pgslw(3);
    pgline(2,[$x_min,$x_max],[$x_min,$x_max]);
    pgsci(1);
    pgsls(1);
    pgslw(1);
}
pgsch(2.5);
$k=0;
#print "$n\n @x\n @y\n";
$nk=0;
for ($i=0;$i<$n;$i++) {
    $X=$x[$i];
    $Y=$y[$i];
    $D[$i]=$X-$Y;
    if ($Y!=0) {
	$R[$k]=$X/$Y;
	$k++;
    }
#    print "$i $X $Y\n";
    pgsci(4);
    pgpoint(1,[$X],[$Y],17);
    pgsci(1);
    pgpoint(1,[$X],[$Y],22);
    if (($X>$x_min)&&($X<$x_max)&&($Y>$y_min)&&($Y<$y_max)) {
	$xx[$nk]=$X;
	$yy[$nk]=$Y;
	$nk++;
    }

}
#print "ONE=$one\n";


#pgpoint($n,\@x,\@y,16);
pgclose;
pgend;


#$med=median(@D);
#$sig=sigma(@D);
#print "DELTA X-Y =$med +- $sig $#D\n";
#$med=median(@R);
#$sig=sigma(@R);
#print "RATIO Y/X =$med +- $sig $#R\n";
$medx=mean(@xx);
$sigx=sigma(@xx);
#print "X =$med +- $sig $nk in the range\n";
$medy=mean(@yy);
$sigy=sigma(@yy);
#print "Y =$med +- $sig $nk in the range\n";
print "$x_min $x_max $y_min $y_max $nk $medx $sigx $medy $sigy\n";

#@REG=regression(\@x,\@y);
#print "REGRESSION @REG\n";

exit;





sub regression {
    local($dx,$dy)=@_;
    my @a=@$dx;
    my @ay=@$dy;
    my $nm=$#a;
    my $reg=Statistics::OLS->new();
    $reg->setData(\@a,\@ay);
    $reg->regress();
    my ($intercept,$slope) = $reg->coefficients();
    my $rsq=$reg->rsq();
    my $r=($rsq)**0.5;
    my ($tstat_intercept,$tstat_slope) = $reg->tstats();
    my $sigma = $reg->sigma();
    my $durbin_watson = $reg->dw();
    if ($rsq==1) {
        $rsq=0.999999;
    }
    my $t = $r*(($n-2)/(1-$rsq))**0.5;
    my $free =$n-2;
    my ($varX,$varY,$covXY) = $reg->var();
    my $stat_x = Math::Stat->new($dx);
    my $average_x = $stat_x->average();
    my $stddev_x = $stat_x->stddev();
    my $stat_y = Math::Stat->new($dy);
    my $average_y = $stat_y->average();
    my $stddev_y = $stat_y->stddev();
    my $S=0;
    my $Sx=0;
    my $Sxx=0;
    my $i=0;
    for ($i=0;$i<$nm;$i++) {
        $Sx=$Sx+($a[$i])**2;
        $Sxx=$Sx+($a[$i])**2;
    }
    my $S=$nm;
    my $delta_a;
    my $delta_b;
    my $delta=$S*$Sxx-($Sx)**2;

    if ($delta!=0) {
        $delta_b = sqrt(abs($S/$delta));
        $delta_a = sqrt(abs($Sxx/$delta));
    }
#    print "$S $nm '$Sx $Sxx' '$delta' '$delta_b' '$delta_a'\n";

    return ($intercept,$delta_b,$slope,$delta_a,$r,$sigma,$t);
}
