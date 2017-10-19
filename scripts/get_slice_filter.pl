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
    print "USE: get_slice_filter.pl INPUT_CUBE.fits filter.dat output.fits [NZ1 NZ2]\n";
    exit;
}

$input_cube=$ARGV[0];
$filter=$ARGV[1];
$out_file=$ARGV[2];



#
#  Reading filter transmission curve
#

$nf=0;
$peak=0;
$wmin=1e12;
$wmax=-1e12;
open(FH,"<$filter");
while ($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if ($#data==2) {
        $flux_f[$nf]=$data[2];
        $wave_f[$nf]=$data[1];
    } else {
        $flux_f[$nf]=$data[1];
        $wave_f[$nf]=$data[0];
    }
    if ($flux_f[$nf]>$peak) {
	$peak=$flux_f[$nf];
    }
    if (($wave_f[$nf]<$wmin)&&($flux_f[$nf]>0)) {
	$wmin=$wave_f[$nf];
    }
    if (($wave_f[$nf]>$wmax)&&($flux_f[$nf]>0)) {
	$wmax=$wave_f[$nf];
    }

    $nf++;
}
close(FH);


$b=rfits("$input_cube");
$h=$b->gethdr;
$crval3=$b->hdr->{CRVAL3};
$cdelt3=$b->hdr->{CDELT3};
$crpix3=$b->hdr->{CRPIX3};
($nx,$ny,$nz)=$b->dims;

if ($#ARGV==4) {
    $NZ1=$ARGV[3];
    $NZ2=$ARGV[4];
} else {
    $NZ1=0;
    $NZ2=$nz-1;
}

#$B=$b->slice(",,$NZ1:$NZ2");
#@N=$B->dims;
#print "$nx,$ny,$nz\n";
#print "@N\n";

for ($i=0;$i<$nz;$i++) {
    $wave[$i]=$crval3+$cdelt3*($i-($crpix3-1));
}


my $pdl_flux = interpol(pdl(@wave), pdl(@wave_f), pdl(@flux_f));
$sum=sum($pdl_flux);
$pdl_flux=$pdl_flux/$sum;

$pdl_out=zeroes($nx,$ny);

for ($i=$NZ1;$i<$NZ2;$i++) {
    $W=$pdl_flux->at($i);
    $a=$b->slice(",,($i)");
    if (($wave[$i]>$wmin)&&($wave[$i]<$wmax)) {
	if (($W>1e-4)) {
	    $pdl_out=$pdl_out+$a*$W;
	}
    }
}
#$pdl_out->sethdr( $h );
for ($i=0;$i<$nx;$i++) {
    for ($j=0;$j<$ny;$j++) {
	$val=$pdl_out->at($i,$j);
	if (abs($val)>1e20) {
	    set($pdl_out,$i,$j,0);
	}
    }
}


$$h{NAXIS}=2;
$pdl_out->inplace->badmask(0);
$pdl_out->sethdr($h);
$pdl_out->wfits($out_file);

exit;






