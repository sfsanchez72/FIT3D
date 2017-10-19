#!/usr/bin/perl
#
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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

if ($#ARGV<2) {
    print "USE: get_slice.pl INPUT_CUBE.fits PREFIX CONF_FILE\n";
    print "CONF_FILE: NAME START_W END_W\n";
    exit;
}

$input_cube=$ARGV[0];
$prefix=$ARGV[1];
$conf_file=$ARGV[2];


$ns=0;
open(FH,"<$conf_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$name[$ns]=$data[0];
	$start_w[$ns]=$data[1];
	$end_w[$ns]=$data[2];
	$ns++;
    }
}
close(FH);

print "$ns slices to cut\n";
print "Reading cube\n";
$b=rfits("$input_cube");
$h=$b->gethdr;
$crval3=$b->hdr->{CRVAL3};
$cdelt3=$b->hdr->{CDELT3};
($nx,$ny,$nz)=$b->dims;

for ($i=0;$i<$ns;$i++) {
    $out_file=$prefix."_".$name[$i]."_".$start_w[$i]."_".$end_w[$i].".fits";
    $start_i=int(($start_w[$i]-$crval3)/$cdelt3);
    $end_i=int(($end_w[$i]-$crval3)/$cdelt3);
    if (($start_i>-1)&&($end_i<$nz)) {
	$npix=$end_i-$start_i+1;
	my $a=$b->slice(",,$start_i:$end_i");
	my $c=average($a->xchg(0,2));
	$c=$c*($end_w[$i]-$start_w[$i]);
	my $d=$c->xchg(0,1);
	$$h{"PIX_WIDTH"}=$npix;
	$$h{"START_W"}=$start_w[$i];
	$$h{"END_W"}=$end_w[$i];
	$$h{"NAXIS"}=2;
#	$c->sethdr( $h );
	$d->wfits($out_file);
	print "$out_file saved\n";
    } else {
	print "section ($start_w[$i],$end_w[$i]) out of margings\n";
    }
}


exit;






