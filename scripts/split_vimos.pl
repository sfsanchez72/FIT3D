#!/usr/bin/perl

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<0) {
    print "USE: split_vimos.pl VIMOS_LR_FITS [Y_SHIFT] \n";
    exit;
}

$infile=$ARGV[0];
$y_shift=0;
if ($#ARGV==1) {
    $y_shift=$ARGV[1];
}

$qua=read_img_header($infile,"HIERARCH ESO DET CHIP1 ID");
if ($qua =~ "59B") {
    $nq=1;
}
if ($qua =~ "59A") {
    $nq=2;
}
if ($qua =~ "60A") {
    $nq=3;
}
if ($qua =~ "60B") {
    $nq=4;
}
print "NQ=$nq\n";


@a_sec0=(800,800,2880,2730);
@b_sec0=(1940,1940,1782,1630);
@c_sec0=(2514,2514,1200,1065);
@d_sec0=(3080,3080,600,460);

#print "NQ=$nq\n";
for ($i=0;$i<4;$i++) {
    $a_sec0[$i]=$a_sec0[$i]+$y_shift;
    $b_sec0[$i]=$b_sec0[$i]+$y_shift;
    $c_sec0[$i]=$c_sec0[$i]+$y_shift;
    $d_sec0[$i]=$d_sec0[$i]+$y_shift;
    $a_sec1[$i]=$a_sec0[$i]+589;
    $b_sec1[$i]=$b_sec0[$i]+589;
    $c_sec1[$i]=$c_sec0[$i]+589;
    $d_sec1[$i]=$d_sec0[$i]+589;
    $a_sec[$i]="".$a_sec0[$i].":".$a_sec1[$i];
    $b_sec[$i]="".$b_sec0[$i].":".$b_sec1[$i];
    $c_sec[$i]="".$c_sec0[$i].":".$c_sec1[$i];
    $d_sec[$i]="".$d_sec0[$i].":".$d_sec1[$i];
}



$pdl=rfits($infile);
$what=":,".$a_sec[$nq-1];
#print "$what\n";
$a_pdl=$pdl->slice($what);
$a_file="a_".$infile;
$a_pdl->hdr->{NAXIS2}=590;
$a_pdl->wfits($a_file);
$what=":,".$b_sec[$nq-1];
#print "$what\n";
$b_pdl=$pdl->slice($what);
$b_file="b_".$infile;
$b_pdl->hdr->{NAXIS2}=590;
$b_pdl->wfits($b_file);
$what=":,".$c_sec[$nq-1];
#print "$what\n";
$c_pdl=$pdl->slice($what);
$c_file="c_".$infile;
$c_pdl->hdr->{NAXIS2}=590;
$c_pdl->wfits($c_file);
$what=":,".$d_sec[$nq-1];
#print "$what\n";
$d_pdl=$pdl->slice($what);
$d_file="d_".$infile;
$d_pdl->hdr->{NAXIS2}=590;
$d_pdl->wfits($d_file);



exit;
