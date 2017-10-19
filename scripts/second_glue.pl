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



$DIR_CONF="/disk-b/sanchez/ppak/legacy/";

open(LOG,">second_glue.log");
open(PROD,">second_glue.out");

$call="rm -f *cube.fits";
mycall($call);
$call="rm -f *sobj.fits";
mycall($call);
$call="rm -f *rss.fits";
mycall($call);
$call="rm -f *.V.fits";
mycall($call);

$call="rm -f SN_*.ps";
mycall($call);
$call="rm -f SN_*.out";
mycall($call);


$call="ls ../pmas*/*.sobj.fits > list.sobj.txt";
mycall($call);

$nobj=0;
$n=0;
open(FH,"<list.sobj.txt");
while($line=<FH>) {
    chop($line);
    @data=split('/',$line);
    $sec=$data[2];
    @data2=split('_',$sec);
    $obj_now=$data2[0];
    $obj[$nobj]=$obj_now;
    $file[$n]=$line;
    $found=0;
    if ($nobj>0) {
	for ($i=0;$i<$nobj;$i++) {
	    if ($obj_now eq $obj[$i]) {
		$found=1;
	    }
	}
    }
    if ($found==0) {
	$nobj++;
    }
    $n++;
}
close(FH);

#print "$n $nobj\n"; exit;

for ($j=0;$j<$nobj;$j++) {
    $n_p_now=0;
    for ($JJ=1;$JJ<4;$JJ++) {
	$check=$obj[$j]."_p".$JJ;
	open(POINT,"cat list.sobj.txt | grep '$check' | wc -l |");
	$point=<POINT>;
	chop($point);
	close(POINT);
	if ($point>0) {
	    $p_now[$n_p_now]=$JJ;
	    $n_p_now++;
	}
    }
    print "$obj[$j] $n_p_now $n\n"; #exit;

    if ($n_p_now==3) {
	for ($i=0;$i<$n;$i++) {
	    @data=split('/',$file[$i]);
	    $sec=$data[2];
	    @data2=split('_',$sec);
	    $obj_now=$data2[0];
	    if ($obj_now eq $obj[$j]) {
		if ($file[$i] =~ "p1") {
		    $call="cp ".$file[$i]." obj_p1.sobj.fits";
		}
		if ($file[$i] =~ "p2") {
		    $call="cp ".$file[$i]." obj_p2.sobj.fits";
		}
		if ($file[$i] =~ "p3") {
		    $call="cp ".$file[$i]." obj_p3.sobj.fits";
		}
#		print "CALL = $call\n"; exit;
		mycall($call);
	    }
	}



	$call="radial_sum_rss.pl obj_p1.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_1.rss.fits";
	mycall($call);
	$call="img2spec.pl rad_1.rss.fits 4 rad_1.txt";
	mycall($call);
	open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_1.txt 0 |");
	while($LINE=<FH>) { 
	    $flux1=$LINE;
	}	chop($flux1);
	chop($flux1);
	close(FH);
	
	$call="radial_sum_rss.pl obj_p2.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_2.rss.fits";
	mycall($call);
	$call="img2spec.pl rad_2.rss.fits 4 rad_2.txt";
	mycall($call);
	open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_2.txt 0 |");
	while($LINE=<FH>) { 
	    $flux2=$LINE;
	}
#	$flux2=<FH>;
	chop($flux2);
	close(FH);
	
	$call="radial_sum_rss.pl obj_p3.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_3.rss.fits";
	mycall($call);
	$call="img2spec.pl rad_3.rss.fits 4 rad_3.txt";
	mycall($call);
	open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_3.txt 0 |");
	while($LINE=<FH>) { 
	    $flux3=$LINE;
	}
#	$flux3=<FH>;
	chop($flux3);
	close(FH);
	
	print "FLUXES = $flux1 $flux2 $flux3\n";
	print LOG "FLUXES = $flux1 $flux2 $flux3\n";
	
	$rat23=$flux2/$flux3;
	if (($rat23<0.95)||($rat23>1.05)) {
	    $call="imarith.pl obj_p2.sobj.fits / ".$rat23." obj_p2.sobj.fits";
	    mycall($call);
	}
	
	$rat13=$flux1/$flux3;
	if (($rat13<0.95)||($rat13>1.05)) {
	    $call="imarith.pl obj_p1.sobj.fits / ".$rat13." obj_p1.sobj.fits";
	    mycall($call);
	}
	
	
#    exit;


	
	$call="Mosaic_rss.pl ".$DIR_CONF."/mos.config mos.rss.fits mos.pt.txt";
	mycall($call);
	$pdl=rfits("mos.rss.fits");
	$cdelt=$pdl->hdr->{CDELT1};
	if ($cdelt==1.5) {
	    $call="rss2rss_mask_ppak.pl mos.rss.fits ".$DIR_CONF."/mos.pt.txt ".$DIR_CONF."/mask.rss.fits mmos.rss.fits";
	} 
	if ($cdelt==3.2) {
	    $call="cp mos.rss.fits mmos.rss.fits";
	}
	print "$call | $cdelt\n";
#	exit;
	mycall($call);
	$call="rm -r ".$obj[$j].".cube.fits ";
	mycall($call);
	$call="rss2cube.tcl mmos.rss.fits mos.pt.txt ".$obj[$j].".cube.fits 3 1e-12 1";
	mycall($call);
	if ($cdelt==3.2) {
	    $call="get_flux_stats_cube.pl ".$obj[$j].".cube.fits 5800 40 3620 3.2 SN_".$obj[$j].".ps/CPS 1 > SN_".$obj[$j].".out";
	}

	if ($cdelt==1.5) {
	    $call="get_flux_stats_cube.pl ".$obj[$j].".cube.fits 5800 40 3845 1.5 SN_".$obj[$j].".ps/CPS 1 > SN_".$obj[$j].".out";
	}

	mycall($call);


	print PROD "$obj[$j].cube.fits $cdelt\n";

	$call="cp mmos.rss.fits mos_".$name_obj[$j].".rss.fits";
	mycall($call);
	$call="cp mos.pt.txt mos_".$name_obj[$j].".pt.txt";
	mycall($call);




	
#    exit;
	
#
# We look for the peak intensity!
#
	$call="get_slice.pl ".$obj[$j].".cube.fits img /disk-b/sanchez/ppak/legacy/slice_V.conf";
	mycall($call);
	$img=rfits("img_V_4500_5500.fits");
	($nx,$ny)=$img->dims;
	$val_max=-1e12;
	$XC=37;
	$YC=37;
	for ($ii=30;$ii<45;$ii++) {
	    for ($jj=30;$jj<45;$jj++) {
		$val=$img->at($ii,$jj);
		if ($val>$val_max) {
		    $val_max=$val;
		    $XC=$ii;
		    $YC=$jj;
		}
	    }
	}
	
	$call="cp img_V_4500_5500.fits ".$obj[$j].".V.fits";
	mycall($call);
	
	$call="rm -r ".$obj[$j].".scube.fits ";
	mycall($call);
	$call="DAR_det_cube ".$obj[$j].".cube.fits  ".$XC." ".$YC." 4 4 ".$obj[$j].".scube.fits 20 1050 3 1 0  0 1";
	mycall($call);

	
	
	print PROD "$obj[$j].scube.fits\n";
	
	
	
    }
    
    
}

$call="rm -f *sobj.fits";
mycall($call);
$call="rm -f rad*rss.fits";
mycall($call);
$call="rm -f mos*";
mycall($call);
$call="rm -f rad*txt";
mycall($call);




exit;

sub mycall {
    my $call=$_[0];
    print "$call | ";
    system($call);
    print LOG "$call\n";
    print "DONE\n";
}

