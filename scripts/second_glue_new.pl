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
$call="rm -f *.B.fits";
mycall($call);

$call="rm -f SN_*.ps";
mycall($call);
$call="rm -f SN_*.out";
mycall($call);


$call="ls ../pmas*/*.sobj.fits > list.sobj.tmp";
mycall($call);


open(FH,"<list.sobj.tmp");
open(OUT,">list.sobj.txt");
while($line=<FH>) {
    chop($line);
    if (($line !~ "vstar")&&($line !~ "obj_")&&($line !~ "run")&&($line !~ "junk")) {
	print OUT "$line\n";
    }
}
close(OUT);
close(FH);


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



print "$n $nobj\n"; 

for ($j=0;$j<$nobj;$j++) {
    print "$j $obj[$j]\n";
}
#    exit;

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



	if (($cdelt==0.7)||($cdelt==2.0)) {
	    $call="radial_sum_rss.pl obj_p1.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 -5 3 rad_1.rss.fits";
	} else {
	    $call="radial_sum_rss.pl obj_p1.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_1.rss.fits";
	}
	mycall($call);
	$call="img2spec.pl rad_1.rss.fits 4 rad_1.txt";
	mycall($call);
	open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_1.txt 0 |");
	while($LINE=<FH>) { 
	    $flux1=$LINE;
	}	chop($flux1);
	chop($flux1);
	close(FH);
	
	if (($cdelt==0.7)||($cdelt==2.0)) {
	    $call="radial_sum_rss.pl obj_p2.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 2 9 rad_2.rss.fits";
	} else {
	    $call="radial_sum_rss.pl obj_p2.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_2.rss.fits";
	}
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
	
    $F=0.05;
    $call="comp_spec.pl rad_2.txt rad_1.txt rad.txt 4000 4500 101";
    mycall($call);
    open(COMP,"<comp_spec.out");
    $comp=<COMP>;
    chop($comp);
    close(COMP);
    @d_comp=split(" ",$comp);
    $rat12=$d_comp[0];
    if (($flux1>0)&&($flux2>0)) {
#	$rat12=$flux2/$flux1;
	if (($rat12<(1-$F))||($rat12>(1+$F))) {
#	    $call="imarith.pl obj_p2.sobj.fits / ".$rat12." obj_p2.sobj.fits";
	    $call="polyfit_spec.pl rad.txt 2 new_rad.txt 300 1900 /null";
	    mycall($call);
	    $call="spec2img.pl new_rad.txt rad.fits 331";
	    mycall($call);
	    $call="imarith.pl obj_p2.sobj.fits '*' rad.fits obj_p2.sobj.fits";
	    mycall($call);
	}
    }
    $call="comp_spec.pl rad_3.txt rad_1.txt rad.txt 4000 4500 101";
    mycall($call);
    open(COMP,"<comp_spec.out");
    $comp=<COMP>;
    chop($comp);
    close(COMP);
    @d_comp=split(" ",$comp);
    $rat13=$d_comp[0];
    if (($flux1>0)&&($flux3>0)) {
#	$rat13=$flux3/$flux1;
	if (($rat13<(1-$F))||($rat13>(1+$F))) {
#	    $call="imarith.pl obj_p3.sobj.fits / ".$rat13." obj_p3.sobj.fits";
	    $call="polyfit_spec.pl rad.txt 2 new_rad.txt 300 1900 /null";
	    mycall($call);
	    $call="spec2img.pl new_rad.txt rad.fits 331";
#	    $call="spec2img.pl rad.txt rad.fits 331";
	    mycall($call);
	    $call="imarith.pl obj_p3.sobj.fits '*' rad.fits obj_p3.sobj.fits";
	    mycall($call);
#	    mycall($call);
	}
    }


#	$rat21=$flux2/$flux1;
#	if (($rat21<0.95)||($rat21>1.05)) {
#	    $call="imarith.pl obj_p2.sobj.fits / ".$rat21." obj_p2.sobj.fits";
#	    mycall($call);
#	}
	
#	$rat31=$flux3/$flux1;
#	if (($rat31<0.95)||($rat31>1.05)) {
#	    $call="imarith.pl obj_p3.sobj.fits / ".$rat31." obj_p3.sobj.fits";
#	    mycall($call);
#	}
	
	
#    exit;


	
	$call="Mosaic_rss.pl ".$DIR_CONF."/mos_new.config mos.rss.fits mos.pt.txt";
	mycall($call);
	$pdl=rfits("mos.rss.fits");
	$cdelt=$pdl->hdr->{CDELT1};
	if ($cdelt==0.7) {
	    $call="rss2rss_mask_ppak.pl mos.rss.fits ".$DIR_CONF."/mos_new.pt.txt ".$DIR_CONF."/mask_V1200.rss.fits mmos.rss.fits 2";
	} 
	if ($cdelt==2) {
	    $call="rss2rss_mask_ppak.pl mos.rss.fits ".$DIR_CONF."/mos_new.pt.txt ".$DIR_CONF."/mask_V500.rss.fits mmos.rss.fits 3";
	} 
	if ($cdelt==1.5) {
	    $call="rss2rss_mask_ppak.pl mos.rss.fits ".$DIR_CONF."/mos.pt.txt ".$DIR_CONF."/mask.rss.fits mmos.rss.fits";
	} 
	if ($cdelt==3.2) {
	    $call="cp mos.rss.fits mmos.rss.fits";
	}
	


	mycall($call);
	$call="cp mmos.rss.fits ".$obj[$j].".rss.fits";
	mycall($call);


#	$call="cp mos.rss.fits mos_".$obj[$j].".rss.fits";
#	mycall($call);
#	$call="cp mos.ps.txt mos_".$obj[$j].".pt.txt";
#	mycall($call);
	
#	exit;


	$call="rm -f ".$obj[$j].".cube.fits ";
	mycall($call);
#	$call="rss2cube.tcl mmos.rss.fits mos.pt.txt ".$obj[$j].".cube.fits 3 1e-12 1";
	$call="rss2cube_int.pl mmos.rss.fits mos.pt.txt 1 ".$obj[$j].".cube.fits 2 5 1";
	mycall($call);

	if ($cdelt==2) {
	    $call="get_flux_stats_cube.pl ".$obj[$j].".cube.fits 5800 40 3745 2.0 SN_".$obj[$j].".ps/CPS 1 > SN_".$obj[$j].".out";
	}
	if ($cdelt==0.7) {
	    $call="get_flux_stats_cube.pl ".$obj[$j].".cube.fits 4500 40 3400 0.7 SN_".$obj[$j].".ps/CPS 1 > SN_".$obj[$j].".out";
	}

	if ($cdelt==3.2) {
	    $call="get_flux_stats_cube.pl ".$obj[$j].".cube.fits 5800 40 3620 3.2 SN_".$obj[$j].".ps/CPS 1 > SN_".$obj[$j].".out";
	}

	if ($cdelt==1.5) {
	    $call="get_flux_stats_cube.pl ".$obj[$j].".cube.fits 5800 40 3845 1.5 SN_".$obj[$j].".ps/CPS 1 > SN_".$obj[$j].".out";
	}

	mycall($call);


	print PROD "$obj[$j].cube.fits $cdelt\n";

	$call="cp mmos.rss.fits mos_".$obj[$j].".rss.fits";
	mycall($call);
	$call="cp mos.pt.txt mos_".$obj[$j].".pt.txt";
	mycall($call);




	
#    exit;
	
#
# We look for the peak intensity!
#

	if ($cdelt!=0.7) {
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
	    $call="DAR_det_cube.pl ".$obj[$j].".cube.fits  ".$XC." ".$YC." 4 4 ".$obj[$j].".scube.fits 20 2050 3 1 0  0 1";
	    mycall($call);

	
	
	    print PROD "$obj[$j].scube.fits\n";
	} else {

	    $call="get_slice.pl ".$obj[$j].".cube.fits img /disk-b/sanchez/ppak/legacy/slice_B.conf";
	    mycall($call);
	    $img=rfits("img_B_3900_4550.fits");
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
	    
	    $call="cp img_B_3900_4550.fits ".$obj[$j].".B.fits";
	    mycall($call);
	    
	    $call="rm -r ".$obj[$j].".scube.fits ";
	    mycall($call);
	    $call="DAR_det_cube.pl ".$obj[$j].".cube.fits  ".$XC." ".$YC." 4 4 ".$obj[$j].".scube.fits 20 2050 3 1 0  0 1";
	    mycall($call);

	}
	
	
    }
    
    
}

$call="rm -f *sobj.fits";
mycall($call);
$call="rm -f rad*rss.fits";
mycall($call);
$call="rm -f mos.*";
mycall($call);
$call="rm -f mmos.*";
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

