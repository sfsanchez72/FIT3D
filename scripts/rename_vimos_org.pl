#!/usr/bin/perl

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

open(DIR,"ls VIMOS.*.fits |");
while($file=<DIR>) {
    chop($file);
    $org=read_img_header($file,"ORIGFILE");
    $call="mv ".$file." ".$org;
    print "$call\n";
    system($call);
}
close(DIR);

exit;
