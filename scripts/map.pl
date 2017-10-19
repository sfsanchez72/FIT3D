#!/usr/bin/perl
#
#
#
use Statistics::OLS;
use Math::Stat;




if ($#ARGV<2) {
    print "USE: map.pl slice.txt/POS_TABLE out.model PREFIX_OUT [OVERRIDE]\n";
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
	if (($n_row==0)&&($over==0)) {
	    print "Low limit for the Flux on system $i:";
	    $cut_flux[$i]=<stdin>;
	    chop($cut_flux[$i]);
	    print "Low limit for the Relative Error of the Flux on system $i:";
	    $cut_error_f[$i]=<stdin>;
	    chop($cut_error_f[$i]);
	    print "Low limit for the Relative Error of the Wavelength on system $i:";
	    $cut_error[$i]=<stdin>;
	    chop($cut_error[$i]);
	    $n_clean[$i]=0;
	}
	$line=<FH>;
	chop($line);
	@data=split(" ",$line);
	$model[$n_row][$i]=$data[0];
	$w[$n_row][$i]=$data[1];
	$e_w[$n_row][$i]=$data[2];
	$f[$n_row][$i]=$data[3];
	$e_f[$n_row][$i]=$data[4];
	$dw[$n_row][$i]=$data[5];
	$e_dw[$n_row][$i]=$data[6];
	if ($model[$n_row][$i] ne "poly1d") {
	    if ($f[$n_row][$i]!=0) {
		if ($over==0) {
		    if (($f[$n_row][$i]>$cut_flux[$i])&&((abs($e_w[$n_row][$i]/$w[$n_row][$i]))<$cut_error[$i])&&((abs($e_f[$n_row][$i]/$f[$n_row][$i]))<$cut_error_f[$i])) {
			$q[$n_row][$i]=1;
			$w_clean[$i][$n_clean[$i]]=$w[$n_row][$i];
			$n_clean[$i]++;
		    } else {
			$q[$n_row][$i]=0;
		    }
		} else {
		    $q[$n_row][$i]=1;
		    $w_clean[$i][$n_clean[$i]]=$w[$n_row][$i];
		    $n_clean[$i]++;
		}
	    } else {
		$q[$n_row][$i]=0;
	    }
	}  else {
		$w_clean[$i][$n_clean[$i]]=0;
		$n_clean[$i]++;
		$q[$n_row][$i]=1;
	}

    }
    $n_row++;

}
close(FH);

for ($i=0;$i<$n_mod;$i++) {
    print "$n_clean[$i] Clean final data for system $i\n";    
    my @x_stat;
    for ($j=0;$j<$n_clean[$i];$j++) {
	$x_stat[$j]=$w_clean[$i][$j];
    }
    $stat = Math::Stat->new(\@x_stat);
    $mean = $stat->average();
    $sig = $stat->stddev();
    print "Sistemic Wavelength of system $i =$mean+-$sig\n";
}

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
	print FH "$k $x[$j] $y[$j] $f[$j][$i] $e_f[$j][$i] $w[$j][$i] $e_w[$j][$i] $dw[$j][$i] $e_dw[$j][$i] $q[$j][$i]\n";
    }
    close(FH);
}

exit;
