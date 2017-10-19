#!/usr/bin/perl

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<0) {
    print "USE: split_vimos_ms_LR.pl VIMOS_LR_FITS \n";
    exit;
}

$infile=$ARGV[0];

$pdl=rfits($infile);
for ($i=0;$i<4;$i++) {
    $k=$i+1;
    if ($k<10) {
	$k="0".$k;
    }
    $name="sec_".$k."_".$infile;
    $start=$i*400;
    $end=($i+1)*400-1;
    $what=":,".$start.":".$end;
    my $a_pdl=$pdl->slice($what);
    $a_pdl->wfits($name);
}
exit;
