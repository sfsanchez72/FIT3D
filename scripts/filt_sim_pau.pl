#!/usr/bin/perl

$w0=3750;
for ($i=0;$i<50;$i++) {
    $name="/home/sanchez/filters/PAU/PAU_".$w0.".txt";
    open(FH,">$name");
    $k=0;
    for ($j=-60;$j<60;$j++) {
	$w=$w0+$j;
	if (($j>-50)&(%j<50)) {
	    $f=1;
	}
	if (abs($j)>=50) {
	    $f=11-0.2*abs($j);
	}
	if (abs($j)>55) {
	    $f=0;
	}
	$k++;
	$f=100*$f;
	print FH "$k $w $f\n";
    }
    close(FH);
    $w0=$w0+100;
}
exit;

