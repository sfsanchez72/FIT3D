#!/usr/bin/perl
#
#
#

#use PGPLOT;

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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");





if ($#ARGV<5) {
    print "USE:fit_spec_back.pl SPECTRUM.TXT CONFIG START_WAVELENGTH END_WAVELENGTH NSIM DEVICE\n";    
    exit;
}

$spec_file=$ARGV[0];
$config=$ARGV[1];
$w_start=$ARGV[2];
$w_end=$ARGV[3];
$nsim=$ARGV[4];
$device=$ARGV[5];

#fit_spec_back.pl sub_spec.txt none 6500 6850 none 0 mask_elines.txt tmp_e.config fit_auto_ssp.UGC9965_cen.out.6500_6850.ps/CPS

$call="fit_spec_back.pl ".$spec_file." none ".$w_start." ".$w_end." none 0 none ".$config." ".$device." > junk.junk";
system($call);
#print "$call\n";
$n=0;
open(FH,"<out.fit_spectra");
while($line=<FH>) {
    chop($line);
    if ($line =~ "eline") {
	@data=split(" ",$line);
        $a_fwhm[$n]=$data[5];
	$n++;

    }
}
close(FH);
$fwhm=mean(@a_fwhm);
$fwhm=3*$fwhm;

open(FH,"<out_mod_res.fit_spectra");
open(MOD,">mod.tmp.txt");
open(RES,">res.tmp.txt");
$n=0;
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$w=$data[0];
	$org=$data[1];
	$mod=$data[2];
	$res=$data[3];
	#print MOD "$n $w $mod\n";
	print MOD "$n $w $org\n";
	print RES "$n $w $mod\n";
	$n++;
    }
}
close(RES);
close(MOD);
close(FH);

$call="montecarlo.pl mod.tmp.txt res.tmp.txt ".$nsim." ".$fwhm." montecarlo.rss.fits";
system($call);
#print "$call\n";

$call="kin_back_rss_auto.pl montecarlo.rss.fits ".$config." none 0  ".$w_start." ".$w_end." fit_spec_back_mc.tmp.out ".$device." > junk.junk";
system($call); #print "$call\n";

#exit;

open(FH,"<fit_spec_back_mc.tmp.out");
for ($i=0;$i<$nsim;$i++) {
    $nmod=<FH>; chop($nmod); 
    for ($j=0;$j<$nmod;$j++) {
	$out=<FH>; chop($out);
        my @data=split(" ",$out);
	$name[$j]=$data[0];
	$nk=0;
	for ($k=0;$k<$#data+1;$k=$k+2) {
	    $flux[$i][$j][$nk]=$data[$k+1];
	    $e_flux[$i][$j][$nk]=$data[$k+2];	   
	    $nk++;
	}
    }
}
close(FH);

#open(OUT,">fit_spec_back_mc.out");
open(OUT,">out.fit_spectra");
print OUT "$nmod\n";
for ($j=0;$j<$nmod;$j++) {
    print OUT "$name[$j] ";
    for ($k=0;$k<$nk;$k++) {
	my @a_flux;
	my @a_e_flux;
	for ($i=0;$i<$nsim;$i++) {
	    $a_flux[$i]=$flux[$i][$j][$k];
	    $a_e_flux[$i]=$e_flux[$i][$j][$k];
	}
	$flux_now=mean(@a_flux);
	$e_flux_now=sigma(@a_flux);
#	print OUT "$flux_now $e_flux_now ";
	printf OUT (" %8.4f %8.4f",$flux_now, $e_flux_now);

    }
    print OUT "\n";
}
close(OUT);



exit;
