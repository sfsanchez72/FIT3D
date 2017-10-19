#!/usr/bin/perl
#
#

if ($#ARGV<3) {
    print "USE: map_arith.pl MAP1 [+/-*] MAP2 MAP3\n";
    exit;
}
$map[0]=$ARGV[0];
$arith=$ARGV[1];
$map[1]=$ARGV[2];
$map[2]=$ARGV[3];



for ($i=0;$i<2;$i++) {
    $n[$i]=0;
    open(FH,"<$map[$i]");
    while ($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$id[$n[$i]][$i]=$data[0];
	$x[$n[$i]][$i]=$data[1];
	$y[$n[$i]][$i]=$data[2];
	$f[$n[$i]][$i]=$data[3];
	$e_f[$n[$i]][$i]=$data[4];
	$w[$n[$i]][$i]=$data[5];
	$e_w[$n[$i]][$i]=$data[6];
	$dw[$n[$i]][$i]=$data[7];
	    $e_dw[$n[$i]][$i]=$data[8];
	$q[$n[$i]][$i]=$data[9];
	$n[$i]++;
    }
    close(FH);
}

if ($n[0]!=$n[1]) {
    print "Files contains different number of data\n";
    exit;
}

open(FHOUT,">$map[2]");
for ($i=0;$i<$n[0];$i++) {
    $q_now=$q[$i][0]*$q[$i][1];
    if ($arith =~ /\-/) {
	$f_now=$f[$i][0]-$f[$i][1];
    }
    if ($arith =~ /\//) {
	if ($f[$i][1]!=0) {
	    $f_now=$f[$i][0]/$f[$i][1];
	} else {
	    $f_now=0.0000;
	    $q_now=0;
	}
    }
    if ($arith =~ /\*/) {
	$f_now=$f[$i][0]*$f[$i][1];
    }
    if ($arith =~ /\+/) {
	$f_now=$f[$i][0]+$f[$i][1];
    }
    $f_now=apr($f_now);
    print FHOUT "$id[$i][0] $x[$i][0] $y[$i][0] $f_now $e_f[$i][0] $w[$i][0] $e_w[$i][0] $dw[$i][0] $e_dw[$i][0] $q_now\n";
}
close(FHOUT);

exit;


sub apr {
    my $z=@_[0];
    my @s=split(/\./,$z);
    my $last=substr($s[1],0,-length($s[1])+4);
    my $zz=$s[0].".".$last;
    return $zz;
}
