#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module
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

use PDL::NiceSlice;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<3) {
    print "USE: rss2rss_mask_ppak.pl INPUT_RSS position_table.txt MASK.fits OUTPUT_RSS [NSHORT]\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$mask_rss=$ARGV[2];
$output_rss=$ARGV[3];
$NS=5;
if ($#ARGV==4) {
    $NS=$ARGV[4];
}


$n=0;
open(PT,"<$pos_table");
$line=<PT>;
chop($line);
@data=split(" ",$line);
$nx_min=1e12;
$ny_min=1e12;
$nx_max=-1e12;
$ny_max=-1e12;

while($line=<PT>) {
    chop($line);
    @data=split(" ",$line);
    $id[$n]=$data[0];
    $xx[$n]=$data[1];
    $yy[$n]=$data[2];
    $n++;
}
close(PT);

$input=rfits($input_rss);
$mask=rfits($mask_rss);
($nx,$ny)=$input->dims;
$output=$input;
$pmask=prodover($mask);
$pmask2=sumover($mask->xchg(0,1));
@data=$pmask2->dims;
#print "@data\n\n";
#exit;


#	$nx1=intval($nx*0.15);
#	$nx2=intval($nx*0.2);
$nx1=$nx;
$val_max=0;
for ($i=0;$i<$nx;$i++) {
    $val=$pmask2->at($i);
#    print "$i $nx1 $val\n";
    if (($i<$nx1)&&($val>($val_max))) {
	$nx1=$i;
	$val_max=$val;
    }
}
$nx2=$nx1+300;
if ($nx2>=$nx) {
    $nx2=$nx-1;
}
if ($nx1==0) {
#    $nx1=intval($nx*0.45);
#    $nx2=intval($nx*0.65);
    $nx1=intval($nx*0.45);
    $nx2=intval($nx*0.55);
}

#print "$nx1,$nx2\n"; exit;

for ($j=0;$j<$ny;$j++) {
    $mask_val=$pmask->at($j);
    if ($mask_val==0) {
	my $t=zeroes($nx);
	my @dist;
	my $nd=0;
	for ($jj=0;$jj<$ny;$jj++) {	    
	    $mask_val=$pmask->at($jj);
	    if ($mask_val!=0) {
		$dist[$jj]=sqrt(($xx[$j]-$xx[$jj])**2+($yy[$j]-$yy[$jj])**2);
	    } else {
		$dist[$jj]=1e12;
	    }

	}

	@dist_sort = sort {$a <=> $b} @dist;
	$ND=0;
	for ($jj=0;$jj<$ny;$jj++) {

	    $mask_val=$pmask->at($jj);
	    if ($dist[$jj]<=$dist_sort[$NS]) {
		$t_now=$input->slice(":,($jj)");
		$t=$t+$t_now;#/(1+$dist[$jj]**2);
		$ND++;
#$ND=$ND+(1+$dist[$jj]**2);
	    }
	
	}

	

	$t=$t/$ND;
#	print "$nx, $nx1,$nx2\n";

	$sum1=medover($input->slice("$nx1:$nx2,($j)"));
	$sum2=medover($t->slice("$nx1:$nx2"));

#	$sum1=sumover($input->slice("$nx1:$nx2,($j)"));
#	$sum2=sumover($t->slice("$nx1:$nx2"));
#	print "$nx, $nx1,$nx2 $sum1,$sum2\n";
	if ($sum2>0) {
	    $rat12=$sum1/$sum2;
	} else {
	    $rat12=1;
	}
	if ($rat12>0) {
	    $t=$t*$rat12;
	}
	for ($i=0;$i<$nx;$i++) {
	    $val=$input->at($i,$j);
	    if ($val==0) {
		$val2=$t->at($i);
		set($output,$i,$j,$val2);
	    }
	}

    $mask_val=$pmask->at($j);
    }


#
# ****************
#
#	set($pmask,$j,1);



#    print "$j/$ny\n";
}
$input->wfits($output_rss);

exit;

sub intval {
    my $a=@_[0];
    my $ia=int($a);
    my $d=$a-$ia;
    if ($d>0.5) {
	$ia=$ia+1;
    }
    return $ia;
}
