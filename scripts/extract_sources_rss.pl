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

if ($#ARGV<4) {
    print "USE: extract_sources_rss.pl INPUT_RSS position_table.txt sources_location.txt OUTPUT_RSS.fits OUTPUT_PT.txt [APERTUIRE]\n";
    print "source_location: ID X Y Aperture\n";
    print "Is APERTURE is defined, its value overrida those in the source_location table\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$loc_table=$ARGV[2];
$output_rss=$ARGV[3];
$out_pos_table=$ARGV[4];
$APER=-1;
if ($#ARGV==5) {
    $APER=$ARGV[5];
}



$n=0;
open(PT,"<$pos_table");
$line=<PT>;
chop($line);
@data=split(" ",$line);
$head=$line;
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
($nx,$ny)=$input->dims;
$h=$input->gethdr;
if ($n!=$ny) {
    print "The number of entries in the position table ($n)\n";
    print "does not correspond with the Y-axis of the RSS ($ny)\n";
    exit;
}


#
# We read the loc table
#

$ns=0;
open(LT,"<$loc_table");
open(OPT,">$out_pos_table");
if ($APER<0) {
    print OPT "$head\n";
} else {
    $A=0.5*$APER;
    print OPT "C $A $A 1\n";
}
while($line=<LT>) {
    chop($line);
    @data=split(" ",$line);
    $ids[$ns]=$data[0];
    $xs[$ns]=$data[1];
    $ys[$ns]=$data[2];
    $aperture[$ns]=$data[3];
    print OPT "$ids[$ns] $xs[$ns] $ys[$ns] 1\n"; 
    if ($APER>$aperture[$ns]) {
	$aperture[$ns]=$APER;
    }
#   print "$ns $aperture[$ns]\n";
    $ns++;
}
close(OPT);
close(LT);
#<stdin>;


$output=zeroes($nx,$ns);
$$h{NAXIS2}=$ns;
$output->sethdr($h);

for ($k=0;$k<$ns;$k++) {
    $first=0;
    for ($i=0;$i<$ny;$i++) {
	$dist=sqrt(($xx[$i]-$xs[$k])**2+($yy[$i]-$ys[$k])**2);
	if ($dist<=$aperture[$k]) {
	    #print "$dist $aperture[$k]\n";
	    $tmp=$output->slice(":,$k");
	    $spec=$input->slice(":,$i");
	    $tmp .=$tmp+$spec;
	    #$tmp2=$output->slice(":,$k");
	    #$tmp2 .= $tmp;

#	    for ($j=0;$j<$nx;$j++) {
#		$val=$input->at($j,$i);
#		$val2=$output->at($j,$k);
#		$tot=$val+$val2;
#		set($output,$j,$k,$tot);
#	    }
	    #if ($first==0) {
	#	$tmp .= $spec;
	#	$first=1;
	#    } else {
#		$tmp .= $tmp+$spec;
	#    }
	}
    }
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
