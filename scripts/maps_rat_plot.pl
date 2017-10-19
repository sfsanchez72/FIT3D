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


if ($#ARGV<4) {
    print "USE: maps_comp_plot.pl map1.fits map2.fits LABEL1 LABEL2 DEV [MIN MAX] [YMIN YMAX] [FIX] [ONE-to-ONE_LINE 0/1]\n";
    exit;
}

$map1=$ARGV[0];
$map2=$ARGV[1];
$label1=$ARGV[2];
$label2=$ARGV[3];
$dev=$ARGV[4];

$x=rfits($map1);
($nx,$ny)=$x->dims;
$n=$nx*$ny;
($mean,$rms,$median,$x_min,$x_max) = stats($x);
print "X=($mean,$rms,$median,$x_min,$x_max)\n";

$y=rfits($map2);
($nx2,$ny2)=$y->dims;
($mean,$rms,$median,$y_min,$y_max) = stats($y);
print "Y=($mean,$rms,$median,$x_min,$x_max)\n";

$y_max=3;
$y_min=-3;


if ($nx2<$nx) {
    $nx=$nx2;
}

if ($ny2<$ny) {
    $ny=$ny2;
}

if ($#ARGV==6) {
    $x_min=$ARGV[5];
    $x_max=$ARGV[6];
}

if ($#ARGV==8) {
    $x_min=$ARGV[5];
    $x_max=$ARGV[6];
    $y_min=$ARGV[7];
    $y_max=$ARGV[8];
}

$fix=0;
if ($#ARGV==9) {
    $x_min=$ARGV[5];
    $x_max=$ARGV[6];
    $y_min=$ARGV[7];
    $y_max=$ARGV[8];
    $fix=$ARGV[9];
}

$one=0;
if ($#ARGV==10) {
    $x_min=$ARGV[5];
    $x_max=$ARGV[6];
    $y_min=$ARGV[7];
    $y_max=$ARGV[8];
    $fix=$ARGV[9];
    $one=$ARGV[10];
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
$kk=0;
#print "$n\n @x\n @y\n";
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$X=$x->at($i,$j);
	
	if ($X>0) {
	    $Y=$y->at($i,$j);
	    $Y=($Y-$X)/$X;
	    $X=log10($X);
	    #$Y=$Y-$X;
	    $XX[$kk]=$X;
	    $YY[$kk]=$Y;
	    $D[$kk]=$X-$Y;
	    $kk++;
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
    }
}
print "ONE=$one\n";


pgclose;
pgend;


$med=median(@D);
$sig=sigma(@D);
print "DELTA X-Y =$med +- $sig $#D\n";

$med=median(@R);
$sig=sigma(@R);
print "RATIO Y/X =$med +- $sig $#R\n";

@REG=regression(\@XX,\@YY);
print "REGRESSION @REG\n";

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
