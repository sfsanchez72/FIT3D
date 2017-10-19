#!/usr/bin/perl


if ($#ARGV<2) {
    print "USE: add_dist.pl dist1.txt dist2.txt out_dist.txt";
    exit;
}

$file1=$ARGV[0];
$file2=$ARGV[1];
$outfile=$ARGV[2];

$n1=0;
open(FH,"<$file1");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id1[$n1]=$data[0];
    $shift1[$n1]=$data[1];
    $n1++;
}
close(FH);

$n2=0;
open(FH,"<$file2");
open(OUT,">$outfile");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id2[$n2]=$data[0];
    $shift2[$n2]=$data[1];
    $shift=$shift1[$n2]+$shift2[$n2];
    print OUT "$id2[$n2] $shift\n";
    $n2++;
}
close(OUT);
close(FH);

exit;
