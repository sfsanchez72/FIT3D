#!/usr/bin/perl

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<0) {
    print "USE: split_vimos_ms_LR.pl VIMOS_LR_FITS\n";
    exit;
}

$infile=$ARGV[0];


for ($i=0;$i<16;$i++) {
    $k=$i+1;
    if ($k<10) {
	$k="0".$k;
    }
    $name="sec_".$k."_".$infile;
    my $a_pdl=rfits($name);
    if ($i==0) {
	($nx,$ny)=$a_pdl->dims;
	$pdl=zeroes($nx,6400);
    }
    $start=$i*400;
    $end=($i+1)*400-1;
    $what=":,".$start.":".$end;
    my $sec=$pdl->slice($what);
    $sec .=$a_pdl;
}
$pdl->hdr->{NAXIS2}=6400;
$pdl->wfits($infile);

exit;
