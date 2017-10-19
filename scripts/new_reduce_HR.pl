#!/usr/bin/perl


if ($#ARGV<7) {
    print "USE: reduce_LR.pl INFILE TRACING_REF WAVELENGTH_REF PLOT N.ADJ EXPTIME EXTINCTION AIRMASS\n";
    exit;
}


@IMG=("B.1","A.2","A.3","B.4");
#@off=(77,-115,-5,54.5,85,-117,-14,53,89,-116.2,-20,55,91,-111.6,0,77.0);
for ($i=0;$i<16;$i++) {
    $off[$i]=0;
}
$off[0]=-1;
$off[3]=-1;

@id=("a","b","c","d");

$infile=$ARGV[0];
$trace=$ARGV[1];
$wave=$ARGV[2];
$plot=$ARGV[3];
$nad=$ARGV[4];
$EXPTIME=$ARGV[5];
$AV=$ARGV[6];
$AIRMASS=$ARGV[7];

$prefix=$infile;
$prefix =~ s/.fits//;
$log="red_".$prefix.".log";

$t_start=time;
#open(FH,"<$config");
open(OUT,">$log");

$call="rm junk?.fits";
action($call);
$call="rm tjunk?.fits";
action($call);
$call="rm wjunk?.fits";
action($call);

for ($i=0;$i<4;$i++) {
    $k=$i+1;
    $call="imarith ".$prefix."_".$IMG[$i].".fits -  gbias".$k.".fits junk".$k.".fits";
    action($call);
    $call="imarith ".$trace."_".$IMG[$i].".fits -  gbias".$k.".fits tjunk".$k.".fits";
    action($call);
    $call="imarith ".$wave."_".$IMG[$i].".fits -  gbias".$k.".fits wjunk".$k.".fits";
    action($call);
}

$n=0;
    for ($j=0;$j<4;$j++) {
	$k=$j+1;
	$OFF=$off[$j];
	$pass="n";
	do {
	    $call="rm tjunk".$k.".trc.fits ";
	    action($call);
	    $call="trace_peaks_recursive tjunk".$k.".fits 1 ".$IMG[$j]."_VIMOS.peaks 10 1 ".$plot." 5 3 tjunk".$k.".trc.fits ".$OFF;
	    action($call);
	    if ($plot==1) {
		print "It is the tracing ok? [y/n]";
		$pass=<stdin>;
		chop($pass);
		if ($pass ne "y") {
		    print "New offset?";
		    $OFF=<stdin>;
		    chop($OFF);
		}
	    } else {
		$pass="y";
	    }

	} while ($pass ne "y");
	$n++;
    }
    


$call="rm junk?.ms.fits";
action($call);
$call="rm junk?.var.fits";
action($call);
$call="rm wjunk?.ms.fits";
action($call);
$call="rm wjunk?.var.fits";
action($call);
$n=0;

    for ($j=0;$j<4;$j++) {
	$k=$j+1;
	if ($nad>0) {
	    $call="extract_gauss_simple junk".$k.".fits 1 tjunk".$k.".trc.fits 7 junk".$k.".ms.fits  junk".$k.".var.fits 3 0 ".$nad;
	} else {
	    $call="extract_aper junk".$k.".fits 1 tjunk".$k.".trc.fits 5 junk".$k.".ms.fits ";
	}
	action($call);
#	if ($nad>0) {
#	    $call="extract_gauss_simple ".$id[$i]."_wjunk".$k.".fits 1 ".$id[$i]."_tjunk".$k.".trc.fits 7 ".$id[$i]."_wjunk".$k.".ms.fits  ".$id[$i]."_wjunk".$k.".var.fits 3 0 ".$nad;
#	} else {
	$call="extract_aper wjunk".$k.".fits 1 tjunk".$k.".trc.fits 5 wjunk".$k.".ms.fits  ";
#	}
	action($call);
    }




$call="glue_vimos_HR_new.pl wjunk ms wjunk.ms.fits";
action($call);


$call="rm wjunk.dc0.fits";
action($call);
$call="dist_cor wjunk.ms.fits  wjunk.dc0.fits wjunk.dist.txt 0 1570 220 2 0";
action($call);
$call="rm wjunk.dc1.fits";
action($call);
$call="rm wjunk.dist.fits";
action($call);
$call="mdist_cor_sp wjunk.dc0.fits 20 50 4 wjunk.dc1.fits wjunk.dist.fits 0 10000 0 30 5 1 ARC_MR.dist.id";
action($call);


$call="rm wjunk.disp_cor.fits";
action($call);
$call="disp_cor wjunk.dc1.fits  4700 2.5 5 50 20 4 wjunk.disp_cor.fits wjunk.disp.txt 0 2865 ARC_MR.disp.id";
action($call);


$call="glue_vimos_HR_new.pl junk ms junk.ms.fits";
action($call);

$call="rm junk.dc1.fits";
action($call);
$call="mdist_cor_external junk.ms.fits wjunk.dist.txt wjunk.dist.fits junk.dc1.fits 0";
action($call);
$call="rm junk.disp_cor.fits";
action($call);
$call="disp_cor_external junk.dc1.fits 4700 2.5 junk.disp_cor.fits wjunk.disp.txt 0";
action($call);


$call="rm junk.FC.fits";
action($call);
$call="imarith junk.disp_cor.fits / fiberflat.fits junk.FC.fits";
#$call="imarith junk.disp_cor.fits / fiberflat.fits ".$prefix.".FC.fits";
action($call);

$call="write_img_header.pl junk.FC.fits CRPIX1 1";
action($call);
$call="write_img_header.pl junk.FC.fits CRVAL1 4700";
action($call);
$call="write_img_header.pl junk.FC.fits CDELT1 2.5";
action($call);

$call="mv junk.FC.fits ".$prefix.".FC.fits";
action($call);

exit;

#$call="rm junk.FC.fits";
#action($call);
#$call="imarith junk.FC0.fits / fiber_trans_cont.fits junk.FC.fits";
#action($call);

if ($EXPTIME>0) {
    $call="flux_calib.pl junk.FC.fits ratio.txt ".$EXPTIME." ".$AV." ".$AIRMASS." junk.fc.fits";
    action($call);
} else {
    $call="mv junk.FC.fits junk.fc.fits";
    action($call);
}



#$call="mv junk.fc.fits ".$prefix.".fc.fits";
#action($call);


#$call="rss2cube_pt_mask.pl junk.fc.fits vimos_pt.txt mask.txt 1 ".$prefix.".cube.fits";
$call="rss2cube_pt_mask.pl junk.fc.fits vimos_pt.txt mask.txt 1 junk.cube.fits";
action($call);

$call="cube2rss.pl junk.cube.fits junk.rss.fits";
action($call);


$call="rm junk.FC.fits";
action($call);
$call="imarith junk.rss.fits / FIBER.rss.fits junk.FC.fits";
action($call);


$call="write_img_header.pl junk.FC.fits CRPIX1 1";
action($call);
$call="write_img_header.pl junk.FC.fits CRVAL1 3700";
action($call);
$call="write_img_header.pl junk.FC.fits CDELT1 5";
action($call);

$call="rss2cube_pt_mask.pl junk.FC.fits VIMOS_CUBE.pt.txt mask.txt 1 ".$prefix.".cube.fits 3";
action($call);

close(OUT);

$delta_t=(time-$t_start)/60;
print "Consumed time: $delta_t\n";
exit;




sub action {
    my $act=$_[0];
    system($act);
    print OUT "$act\n";
    print "$act\n";
    return;
}



$n=0;
for ($i=0;$i<4;$i++) {
    for ($j=0;$j<4;$j++) {
	$k=$j+1;
	$call="dist_cor ".$id[$i]."_wjunk".$k.".ms.fits ".$id[$i]."_wjunk".$k.".dc0.fits ".$id[$i]."_wjunk".$k.".dist.txt 0 420 60 2 0 30 100 400";
	action($call);
	$call="mdist_cor_sp ".$id[$i]."_wjunk".$k.".dc0.fits 3 40 4 ".$id[$i]."_wjunk".$k.".dc1.fits ".$id[$i]."_wjunk".$k.".disp.fits  0 3000 1 30 5 10";
	action($call);
    }
}

