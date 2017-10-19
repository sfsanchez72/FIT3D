#!/usr/bin/perl

$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");


open(DIR,"ls VIMOS*.fits |");
open(OUT,">quads.txt");
while($file=<DIR>) {
    chop($file);
    $qua=read_img_header($file,"HIERARCH ESO DET CHIP1 ID");
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
    print "$file $nq\n";
    print OUT "$file $nq\n";
}
close(OUT);
close(DIR);

exit;
