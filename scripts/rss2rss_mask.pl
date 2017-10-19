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
    print "USE: rss2rss_map.pl INPUT_RSS position_table.txt MASK.fits OUTPUT_RSS\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$mask_rss=$ARGV[2];
$output_rss=$ARGV[3];

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

$dlimit=6;
for ($i=0;$i<$nx;$i++) {
    $sec=$mask->slice("($i),");
    $p=prodover($sec);
    $s=sumover($sec);
    if (($s>0)&&($p==0)) {
	for ($j=0;$j<$ny;$j++) {
	    $val_mask=$mask->at($i,$j);
	    if ($val_mask==0) {
		my @out;
		$nout=0;
		$sum=0;
		$prod=1;
		$dmin=1e12;
		for ($jj=0;$jj<$ny;$jj++) {
		    $dist=sqrt(($xx[$j]-$xx[$jj])**2+($yy[$j]-$yy[$jj])**2);
		    $mask_now=$mask->at($i,$jj);
		    if (($dist<$dlimit)&&($mask_now==1)) {
			if ($dist<$dmin) {
			    $dmin=$dist;
			}
			$val=$input->at($i,$jj);
			$out[$nout]=$val;
			$nout++;
		    }
		}	    
		if ($nout>0) {
		    $val=$input->at($i,$j);
		    $med=median(@out);
#		    print "$nout $med $val $i,$j $jj $dist $dlimit $dmin ($nx,$ny)\n";
		    
		    set($output,$i,$j,$med);
		}
	    }
	}
    }
    print "$i/$nx\n";
}

$output->wfits($output_rss);

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
