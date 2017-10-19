#!/usr/bin/perl
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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";


if ($#ARGV<6) {
    print "USE read_ised_ASCII.pl ised_ASCII CRVAL CDELT NPIX AGE_index(0/221) output_spec.txt MEDIAN_BOX\n";
    exit;
}

$infile=$ARGV[0];
$crval=$ARGV[1];
$cdelt=$ARGV[2];
$npix=$ARGV[3];
$age_index=$ARGV[4];
$outfile=$ARGV[5];
$med_box=$ARGV[6];

$nsec=1221;
if ($infile =~ "_hr_") {
    $nsec=6900;
}

#print "NSEC=$nsec\n";

my @input;
$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    for ($i=0;$i<($#data+1);$i++) {
	$input[$n]=$data[$i];
#	print "$n $input[$n]\n";
	$n++;
    }
#    $input=$input." ".$line;
}
#print "DONE\n";
close(FH);

$nsteps=$input[0];
#print "nsteps=$nsteps\n";
for ($i=0;$i<$nsteps;$i++) {
    $tb[$i]=$input[1+$i]/1e9;
#    print "$i $tb[$i]\n";
}
$ml=$input[1+$nsteps];
$mo=$input[2+$nsteps];
$iseg=$input[3+$nsteps];
#print "ml=$ml, mo=$mo, iseg=$iseg\n";
$ns=4+$nsteps;
for ($i=0;$i<$iseg;$i++) {
    $xx[$i]=$input[$ns];
    $lm[$i]=$input[$ns+1];
    $um[$i]=$input[$ns+2];
    $baux[$i]=$input[$ns+3];
    $cn[$i]=$input[$ns+4];
    $cc[$i]=$input[$ns+5];
#    print "$i $xx[$i] $lm[$i] $um[$i] $baux[$i] $cn[$i] $cc[$i]\n";
    $ns=$ns+6;
}
#$ns=$ns
$totm=$input[$ns];
$totn=$input[$ns+1];
$avs=$input[$ns+2];
$tau=$input[$ns+3];
$tau1=$input[$ns+4];
$tau2=$input[$ns+5];
$tau3=$input[$ns+6];
print "$totm $tau $tau1 $tau2 $tau3\n";
$i=0;
while ($input[$i]!=$nsec) {
    $i++;
}

for ($nm=0;$nm<$nsteps;$nm++) {
    $i++;
    for ($ni=0;$ni<$nsec;$ni++) {
	$flux[$nm][$ni]=$input[$i];
	$i++;
    }
    while ($input[$i]!=$nsec) {
	$i++;
    }
}

for ($j=0;$j<$ni;$j++) {
    $wave_in[$j]=$flux[0][$j];
    $flux_in[$j]=$flux[$age_index][$j];
}

for ($j=0;$j<$npix;$j++) {
    $wave_out[$j]=$crval+$cdelt*$j;    
}


my $flux_out = interpol(pdl(@wave_out), pdl(@wave_in), pdl(@flux_in));

if ($med_box>1) {
    $a=ones($med_box);
    $flux_out=$flux_out->conv1d($a);
}


open(OUT,">$outfile");
for ($j=0;$j<$npix;$j++) {
    $k=$j+1;
    $val=$flux_out->at($j);
#    print OUT "$k $flux[0][$j] $flux[$age_index][$j]\n";
    print OUT "$k $wave_out[$j] $val\n";
}
close(OUT);

exit;
