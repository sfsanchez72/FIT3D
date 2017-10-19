#!/usr/bin/perl
#
#
#
use Statistics::OLS;
use Math::Stat;




if ($#ARGV<2) {
    print "USE: map_back.pl slice.txt out.model PREFIX_OUT [OVERRIDE]\n";
    exit;
}

$slice=$ARGV[0];
$out_model=$ARGV[1];
$prefix=$ARGV[2];
$over=0;
if ($#ARGV==4) {
    $over=$ARGV[3];
} 


$ns=0;
open(FH,"<$slice");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
#    print "$x[$ns] $y[$ns]\n";
	$ns++;
    }
}
close(FH);

$n_row=0;
open(FH,"<$out_model");
while($line=<FH>) {
#    chop($line);
#    $n_row=$line;
#    $line=<FH>;
    chop($line);
    $n_mod=$line;
    for ($i=0;$i<$n_mod;$i++) {
	$line=<FH>;
	chop($line);
	@data=split(" ",$line);
	$model[$n_row][$i]=$data[0];
	$w[$n_row][$i]=$data[1];
	$e_w[$n_row][$i]=$data[2];
	$f[$n_row][$i]=$data[3];
	$e_f[$n_row][$i]=$data[4];
	$dw[$n_row][$i]=2.345*$data[5];
	$e_dw[$n_row][$i]=2.345*$data[6];
	$v[$n_row][$i]=$data[7];
	$e_v[$n_row][$i]=$data[8];
	$AV[$n_row][$i]=$data[9];
	$e_AV[$n_row][$i]=$data[10];

    }
    $n_row++;

}
close(FH);

for ($i=0;$i<$n_mod;$i++) {
    if($i<10) {
	$num="0".$i;
    } else {
	$num=$i;
    }
    $outfile=$prefix."_".$num.".map";
    open(FH,">$outfile");
    for ($j=0;$j<$n_row;$j++) {
	$k=$j+1;
	print FH "$k $x[$j] $y[$j] $f[$j][$i] $e_f[$j][$i] $v[$j][$i] $e_v[$j][$i] $dw[$j][$i] $e_dw[$j][$i] $w[$j][$i] $e_w[$j][$i] $q[$j][$i] $AV[$j][$i] $eAV[$j][$i]\n";
    }
    close(FH);
}

exit;
