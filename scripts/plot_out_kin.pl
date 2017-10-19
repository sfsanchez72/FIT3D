#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
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

$galfit="/work1/ssa/galfit/galfit";
$cl="/work1/ssa/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<3) {
    print "USE: plot_out_kin.pl INPUT.org INPUT.mod INPUT.res DEVICE [pos_table.txt] [MIN MAX]\n";
    exit;
}

$input_org=$ARGV[0];
$input_mod=$ARGV[1];
$input_res=$ARGV[2];
$dev=$ARGV[3];
$tab=0;
if ($#ARGV==4) {
    $pos_table=$ARGV[4];
    $tab=1;
}

$def=0;
if ($#ARGV==6) {
    $min=$ARGV[5];
    $max=$ARGV[6];
    $def=1;
}



$org=rfits($input_org);

if ($tab==1) {
    open(FH,"<$pos_table");
    open(OUT,">slice.tmp");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	    $id[$ns]=$data[0];
	    $x[$ns]=$data[1];
	    $y[$ns]=$data[2];
	    $slice=$org->slice(",$ns");
	    $flux[$ns]=sum($slice);
	    $k=$ns+1;
	    print OUT "$k $x[$ns] $y[$ns] $flux[$ns]\n";
	    $ns++;
	} else {
	    $shape=$data[0];
	    $size=$data[1]*2;
	}
    }
    close(OUT);
    close(FH);
}



$mod=rfits($input_mod);
$res=rfits($input_res);

($nx,$ny)=$org->dims();
$crval=$org->hdr->{CRVAL1};
$cdelt=$org->hdr->{CDELT1};
$cpix=$org->hdr->{CRPIX1};

$command="";
$n=0;
while ($command ne "X") {
    print "Plotted Spectrum Number $n\n";
    print ") => Next Spectra\n";
    print "( => Previous Spectra\n";
    print "n => Number of spectrum [0,$ny]\n";
    print "X => Exit\n";
    open(OUT,">plot_out_kin.spec");
    for ($i=0;$i<$nx;$i++) {
	$w=$crval+$cdelt*($i+($cpix-1));
	$org_val=$org->at($i,$n);
	$res_val=$res->at($i,$n);
	$mod_val=$mod->at($i,$n);
	print OUT "$w $org_val $mod_val $res_val 0\n";
    }
    close(OUT);
    $call="plot_out_fit.pl plot_out_kin.spec 1".$dev." ";
    if ($def==1) {
	$call=$call." ".$min." ".$max;
    }
    system($call);
#    $command=<stdin>;
#    chop($command);
#    $command=getc(stdin);
#    $command="";
#    while ($command eq "") {
#    while (($command=getc) ) {
#	sysread(STDIN,$command,1);
#    }
    if ($tab==1) {
	$call="plot_slice.pl slice.tmp 2".$dev." ".$shape." ".$size." ".$n;
	system($call);
    }

    $command=getc;
    if ($command eq "(") {
	$n=$n-1;
	if ($n<0) {
	    $n=$ny-1;
	}
    }
    if ($command eq "\n") {
	$n=$n+1;
	if ($n>=($ny)) {
	    $n=0;
	}
    }
    if ($command eq ")") {
	$n=$n+1;
	if ($n>=($ny)) {
	    $n=0;
	}
    }
    if ($command eq "n") {	    
	$val=<stdin>;
	chop($val);
	$n=floor($val);
	if ($n>=($ny)) {
	    $n=0;
	}
	if ($n<0) {
	    $n=$ny-1;
	}
    }

}

exit;    ;
