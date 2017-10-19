#!/usr/bin/perl
#
#
#
use PGPLOT;


use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");

# This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s the and output file
# like the one created by "fit_3D_fg.pl"
# and fix one parameter to do
# a force fitting latter
if ($#ARGV<2) {
    print "USE: compare_param.pl FILE_PARAM1 FILE_PARAM2 N_SIGMA\n";        
    exit;
}

$file1=$ARGV[0];
$file2=$ARGV[1];
$nsigma=$ARGV[2];

#
# Reading the input file
#

$n1=0;
open(INPUT,"<$file1");
while($input=<INPUT>) {
    chop($input);
    @data=split(" ",$input);
    $id1[$n1]=$data[0];
    $flux1[$n1]=$data[1];
    $n1++;
};
close(INPUT);

$n2=0;
open(INPUT,"<$file2");
while($input=<INPUT>) {
    chop($input);
    @data=split(" ",$input);
    $id2[$n2]=$data[0];
    $flux2[$n2]=$data[1];
    $n2++;
};
close(INPUT);

if ($n1!=$n2) {
    print "WARNING $n1!=$n2\n";
}

for ($i=0;$i<$n1;$i++) {
    $dif[$i]=$flux2[$i]-$flux1[$i];
}
$mean=mean(@dif);
$sigma=sigma(@dif);
$median=median(@dif);
print "$mean ($median) +- $sigma ($n1)\n";

$j=0;
for ($i=0;$i<$n1;$i++) {
    if (abs($dif[$i]-$mean)<($nsigma*$sigma)) { 
	$dif2[$j]=$flux2[$i]-$flux1[$i];
	$j++;
    }
}
$mean2=mean(@dif2);
$sigma2=sigma(@dif2);
$median2=median(@dif2);
print "$mean2 ($median2) +- $sigma2 ($j)\n";
exit;
