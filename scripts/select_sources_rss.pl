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

if ($#ARGV<5) {
    print "USE: select_sources_rss.pl INPUT_RSS position_table.txt MIN_WAVELENGTH MAX_WAVELENGTH NSIGGIMA sources_location.txt [CONTINUM_SUB=0/1]\n";
    print "source_location: ID X Y Aperture\n";
#    print "Is APERTURE is defined, its value overrida those in the source_location table\n";
    exit;
}

$input_rss=$ARGV[0];
$pos_table=$ARGV[1];
$min_wave=$ARGV[2];
$max_wave=$ARGV[3];
$nsigma=$ARGV[4];
$loc_table=$ARGV[5];
$c_sub=0;
if ($#ARGV==6) {
    $c_sub=1;
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

#print "N=$n\n";

$input=rfits($input_rss);
($nx,$ny)=$input->dims;
if ($n!=$ny) {
    print "The number of entries in the position table ($n)\n";
    print "does not correspond with the Y-axis of the RSS ($ny)\n";
    exit;
}


$h=$input->gethdr;
$crval=$$h{CRVAL1};
$cdelt=$$h{CDELT1};
$crpix=$$h{CRPIX1};
$n1=int(($min_wave-$crval)/$cdelt)+$crpix-1;
$n2=int(($max_wave-$crval)/$cdelt)+$crpix-1;
$n0=$n1-20;
$n3=$n2+20;
if ($n0<0) {
    $n0=0;
}
if ($n2>=$nx) {
    $n2=$nx-1;
}

my $a=$input->slice("$n1:$n2,");
my $c=average($a);
if ($c_sub==1) {
    my $aC1=$input->slice("$n0:$n1,");
    my $cC1=average($aC1);
    my $aC2=$input->slice("$n2:$n3,");
    my $cC2=average($aC2);
    $aC=0.5*($cC1+$cC2);
    $c=$c-$aC;
}


$c=$c*($max_wave-$min_wave);

#print "$c\n";

$med=median(list($c));
$sig=sigma(list($c));

for ($niter=0;$niter<3;$niter++) {
    my @out;
    my $nout=0;
    for ($i=0;$i<$n;$i++) {
#	print "$niter $i\n";
	$val=$c->at($i);
#	print "$val **\n";
	if ($val<($med+$nsigma*$sig)) {
	    $out[$nout]=$val;
	    $nout++;
	}	
    }
    $med=median(@out);
    $sig=sigma(@out);
}




open(FH,">$loc_table");
for ($i=0;$i<$n;$i++) {
    $val=$c->at($i);
    if ($val>($med+$nsigma*$sig)) {
	print FH "$i $xx[$i] $yy[$i] 3\n";
    }	
}
close(FH);

exit;


#
# We read the loc table
#

$ns=0;
open(LT,"<$loc_table");
open(OPT,">$out_pos_table");
print OPT "$head\n";
while($line=<LT>) {
    chop($line);
    @data=split(" ",$line);
    $ids[$ns]=$data[0];
    $xs[$ns]=$data[1];
    $ys[$ns]=$data[2];
    $aperture[$ns]=$data[2];
    print OPT "$ids[$ns] $xs[$ns] $ys[$ns] 1\n"; 
    if ($APER>$aperture[$ns]) {
	$aperture[$ns]=$APER;
    }
    $ns++;
}
close(OPT);
close(LT);

$output=zeroes($nx,$ns);
$$h{NAXIS2}=$ns;
$output->sethdr($h);

for ($k=0;$k<$ns;$k++) {
    my $tmp=$output->slice(":,$k");
    for ($i=0;$i<$ny;$i++) {
	$dist=sqrt(($xx[$i]-$xs[$k])**2+($yy[$i]-$ys[$k])**2);
	if ($dist<=$aperture[$k]) {
	    my $spec=$input->slice(":,$i");
	    $tmp .= $tmp+$spec;
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
