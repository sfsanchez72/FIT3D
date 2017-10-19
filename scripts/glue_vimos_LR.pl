#!/usr/bin/perl

use PDL::Core;

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<2) {
    print "USE: glue_vimos_LR.pl VIMOS_LR_FITS[1-4].fits type(e.g., ms) OUT_FILE.fits\n";
    exit;
}

$sufix=$ARGV[0];
$type=$ARGV[1];
$outfile=$ARGV[2];


@ord=('d','c','b','a');
$n_img=0;
for ($n=4;$n>0;$n--) {
    for ($i=0;$i<4;$i++) {
	$infile=$ord[$i]."_".$sufix."".$n.".".$type.".fits";
	my $in_pdl=rfits($infile);
	if (($n==4)&&($i==0)) {
	    $hdr=$in_pdl->gethdr();
	    $naxis1=$hdr->{NAXIS1};
	    $pdl=zeroes($naxis1,6400);
	    $pdl->sethdr($hdr);
	    $pdl->hdr->{NAXIS2}=6400;
	}
	$j_start=$n_img*400;
	$j_end=$j_start+400;
	$jj=0;
	for ($j=$j_start;$j<$j_end;$j++) {
	    for ($k=0;$k<$naxis1;$k++) {
#		print "($k,$j) ($k,$jj)\n";
		set($pdl,$k,$j,$in_pdl->at($k,$jj));
	    }
	    $jj++;
	}
	print "$infile included\n";
	$n_img++;
    }
}

$pdl->wfits($outfile);
exit;
