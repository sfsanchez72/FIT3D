#!/usr/bin/perl
#
#
#

print "-----------------------------------\n";
print "@ S.F.Sanchez, August 2006, CAHA\n";
print "-----------------------------------\n";
print "This program is thought to perform a quick\n";
print "reduction of PMAS data to help you on visualizing the \n";
print "results using R3D. It is not thought to produce final \n";
print "products, and we strongly recomend to use it at the\n";
print "telescope for a visual inspection of the data\n";
print "It does only perform a 0-order correction of the distorsion\n";
print "but it can help you to visualize the data\n";
print "-----------------------------------\n";

print "********************\n";
print "Are you using PPAK(0) or the LArr(1) mode?\n";
$type=<stdin>;
chop($type);
if ($type==0) {
    $n_exp=382;
    $n_fib=2;
} else {
    $n_exp=256;
    $n_fib=0;
}
print "N.FIBERS=$n_exp\n";

print "********************\n";
print "Input RAW frame to be reduced?\n";
$infile=<stdin>;
chop($infile);
($prefix,$junk)=split(/\./,$infile);


print "********************\n";
print "BIAS frame [ENTER for no bias subtraction]?\n";
$bias=<stdin>;
chop($bias);

if ($bias ne "") {
    $call="rm junk.fits";
    print "$call\n";
    system($call);
    $call="imarith ".$infile." - ".$bias." junk.fits";
    print "$call\n";
    system($call);
} else {
    $call="cp ".$infile." junk.fits";
    print "$call\n";
    system($call);
}

print "********************\n";
print "CONTINUUM frame to find & trace the spectra?\n";
$trace=<stdin>;
chop($trace);

if ($trace ne "") {
    $call="rm tjunk.fits";
    print "$call\n";
    system($call);
    $call="imarith ".$trace." - ".$bias." tjunk.fits";
    print "$call\n";
    system($call);
    $n_found=0;
    $limit=0.05;
#    $limit=0.5;
    do {
	$call="peak_find tjunk.fits 0 1 1 5 1 3 ".$limit." tjunk.peaks";
	print "$call\n";
	system($call);
	print "NPEAKS_EXPECTED=$n_exp\n";
	open(FH,"wc -l tjunk.peaks |");
	$line=<FH>;
	@data=split(" ",$line);
	$n_found=$data[0]-1;
	close(FH);
	if ($n_found!=$n_exp) {
	    print "********************\n";
	    print "We need to redefine the thershold (green-line)\n";
	    print "New limit (0.01-1.00)?";
	    $limit=<stdin>;
	    chop($limit);
	}
    } while ($n_found!=$n_exp);
    $call="rm tjunk.trc.fits";
    print "$call\n";
    system($call);
    $call="trace_peaks_recursive tjunk.fits 0 tjunk.peaks 1 3 0 5 2 tjunk.trc.fits 0 2.5";
    print "$call\n";
    system($call);
    $call="rm junk.ms.fits";
    print "$call\n";
    system($call);
    $call="extract_aper junk.fits 0 tjunk.trc.fits 5 junk.ms.fits"; 
    print "$call\n";
    system($call);
} else {
    print "We cannot continue the reduction\n";
    exit;
}


print "********************\n";
print "ARC-LAMP frame to correct for distorsion & dispersion?\n [ENTER if not correction $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #d]\n";
$arc=<stdin>;
chop($arc);

if ($arc ne "") {
    $call="rm arc.fits";
    print "$call\n";
    system($call);
    $call="imarith ".$arc." - ".$bias." arc.fits";
    print "$call\n";
    system($call);
    $call="rm arc.ms.fits";
    print "$call\n";
    system($call);
    $call="extract_aper arc.fits 0 tjunk.trc.fits 5 arc.ms.fits"; 
    print "$call\n";
    system($call);

    print "********************\n";
    print "You need to select the pixel of an emission line and a certain width\n";
    print "to trace the emission line along the cross-dispersion or pseudo-list axis (X-axis)\n";
    $end=0;
    while ($end==0) {
	$call="image_plot.pl arc.ms.fits /xs"; 
	print "$call\n";
	system($call);
	print "Pixel of an emission line?";
	$start=<stdin>;
	chop($start);
	print "Width of the emission line search?";
	$width=<stdin>;
	chop($width);
	$call="rm arc.dc0.fits";
	print "$call\n";
	system($call);    
	$call="dist_cor arc.ms.fits arc.dc0.fits arc.dist.txt 0 ".$start." ".$width." 2 1 ".$n_fib;
	print "$call\n";
	system($call);
	$call="image_plot.pl arc.dc0.fits /xs"; 
	print "$call\n";
	system($call);    
	print "Are you happy with the 0-order distorsion correction? [0 or press=NO, 1=YES]";
	$end=<stdin>;
	chop($end);
    }
    $call="rm junk.dc0.fits";
    print "$call\n";
    system($call);    
    $call="dist_cor_external junk.ms.fits arc.dist.txt junk.dc0.fits 0";
    print "$call\n";
    system($call);    
    $call="cp junk.dc0.fits ".$prefix.".disp_cor.fits";
    print "$call\n";
    system($call);    

    print "Do you want to perform a second order correction [1=YES,0 or ENTER=NO]?";
    $second=<stdin>;
    chop($second);

    if ($second==1) {
	$end=0;
	while ($end==0) {
	    $call="rm arc.dc1.fits";
	    print "$call\n";
	    system($call);    
	    $call="rm arc.dist.fits";
	    print "$call\n";
	    system($call);    
	    if ($type==0) {
		$call="mdist_cor_sp arc.dc0.fits 5 30 4 arc.dc1.fits arc.dist.fits 22 10100 1 0 0 2";
		print "$call\n";
		system($call);    
	    } else {
		$call="mdist_cor_sp arc.dc0.fits 5 30 4 arc.dc1.fits arc.dist.fits 22 10100 1";
		print "$call\n";
		system($call);    
	    }
	    $call="image_plot.pl arc.dc1.fits /xs"; 
	    print "$call\n";
	    system($call);    
	    print "Are you happy with the 1-order distorsion correction? [0 or press=NO, 1=YES]";
	    $end=<stdin>;
	    chop($end);
	}
	$call="rm junk.dc1.fits";
	print "$call\n";
	system($call);    
	$call="mdist_cor_external junk.ms.fits arc.dist.txt arc.dist.fits junk.dc1.fits 0";
	print "$call\n";
	system($call);    
	$call="cp junk.dc1.fits ".$prefix.".disp_cor.fits";
	print "$call\n";
	system($call);   
    }
    
#mdist_cor_sp glist_all.dc0.fits 2 30 4 glist_all.dc1.fits glist_all.dist.fits 22 10100 1 0 0 2

    
} else {
    $call="cp junk.ms.fits ".$prefix.".disp_cor.fits";
    print "$call\n";
    system($call);    
}

$call="write_img_header.pl ".$prefix.".disp_cor.fits CRPIX1 1";
print "$call\n";
system($call);    
print "CRVAL?\n";
$crval=<stdin>;
chop($crval);
$call="write_img_header.pl ".$prefix.".disp_cor.fits CRVAL1 ".$crval;
print "$call\n";
system($call);    
print "CDELT?\n";
$cdelt=<stdin>;
chop($cdelt);
if ($cdelt<=0) {
    $cdelt=1;
}
$call="write_img_header.pl ".$prefix.".disp_cor.fits CDELT1 ".$cdelt;
print "$call\n";
system($call);    

$call="ps -fea | grep tk_e3d > e3d.active";
system($call);

open(FH,"<e3d.active");
$n=0;
while($line=<FH>) {
    $n++;
}
close(FH);
#print "NN=$n\n";
if ($n==2) {
    $call="tk_e3d.tcl &";
    print "$call\n";
    system($call);    
    print "Press Enter to continue";
    <stdin>;
}

open(FH,"pwd |");
$pwd=<FH>;
chop($pwd);
close(FH);

open(FH,"which tk_e3d |");
$e3d_dir=<FH>;
chop($e3d_dir);
close(FH);
$clean="user/bin/tk_e3d";
$e3d_dir =~ s/$clean//;

if ($type==0) {
    $call="split_ppak.pl ".$prefix.".disp_cor.fits ".$prefix.".obj.fits cal.fits sky.fits 0";
    print "$call\n";
    system($call);    
    $call="write_img_header.pl ".$prefix.".obj.fits CRVAL1 ".$crval;
    print "$call\n";
    system($call);    
    $call="write_img_header.pl ".$prefix.".obj.fits CDELT1 ".$cdelt;
    print "$call\n";
    system($call);    
    $e3d="new_import_rss ".$pwd."/".$prefix.".obj.fits ".$e3d_dir."/data/ppak_pt_arc.txt";

} else {
    $e3d="new_import_rss ".$pwd."/".$prefix.".disp_cor.fits ".$e3d_dir."/data/pmas_pt.txt";
}
#print "E3D=$e3d\n";

open(FH,">input.e3d");
print FH "$e3d\n";
print FH "euro3d ask_e3d_info\n"; 
print FH "euro3d recolor_image grey 1 0.65 0.65 1.0 \n";
print FH "euro3d draw_raw 22.375869 106.387664 1 0 grey 0.65 0.65 1.0 \n";
close(FH);

system("cp input.e3d ~/.E3D/");



exit;


