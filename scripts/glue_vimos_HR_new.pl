#!/usr/bin/perl

use PDL::Core;

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<2) {
    print "USE: glue_vimos_HR.pl VIMOS_LR_FITS[1-4].fits type(e.g., ms) OUT_FILE.fits\n";
    exit;
}

$sufix=$ARGV[0];
$type=$ARGV[1];
$outfile=$ARGV[2];


#@ord=('d','c','b','a');
$n_img=0;
for ($n=4;$n>0;$n--) {
#    for ($i=0;$i<4;$i++) {
#	$infile=$ord[$i]."_".$sufix."".$n.".".$type.".fits";
	$infile=$sufix."".$n.".".$type.".fits";
#	print "INFILE=$infile\n";
	my $in_pdl=rfits($infile);
	if ($n==4) {
	    $hdr=$in_pdl->gethdr();
	    $naxis1=$hdr->{NAXIS1};
	    $pdl=zeroes($naxis1,1600);
	    $pdl->sethdr($hdr);
	    $pdl->hdr->{NAXIS2}=1600;
	}
	$j_start=$n_img*400;
	$j_end=$j_start+400-1;

 #	$jj=0;
#	for ($j=$j_start;$j<$j_end;$j++) {
#	    for ($k=0;$k<$naxis1;$k++) {
#		set($pdl,$k,$j,$in_pdl->at($k,$jj));
#	    }
#	    $jj++;end
#	}
	if ($n!=2) {
	$t=$pdl->slice(":,$j_start:$j_end"); 
	$t .=$in_pdl;#->slice(":,0:400"));
	} else {
		$jjj=0;
	for ($j=$j_start;$j<$j_end;$j=$j+20) {
		$jj=$j+19;
		$t=$pdl->slice(":,$jj:$j"); 
		$jjjj=$jjj+19;
		$t .=$in_pdl->slice(":,$jjj:$jjjj");
		$jjj=$jjj+20;
		}
	}


	print "$infile included $n_img [$j_start,$j_end]\n";
	$n_img++;
#    }
}

$pdl->wfits($outfile);
exit;
