#!/usr/bin/perl
#
#




open(FH,">pt_com.txt");
print FH "R 0.5 0.5 0\n";
$n=1;
for ($j=0;$j<36;$j++) {
    for ($i=0;$i<33;$i++) {
	$x=0.5*$j;
	$y=0.5*$i;
	$k=36-$j;
	$type=1;
	if (($k>16)&&($i<8)) {
	    $type=0;
	}
	if (($k>26)&&($i<17)) {
	    $type=0;
	}
	if (($k<11)&&($i>15)) {
	    $type=0;
	}
	if (($k<21)&&($i>23)) {
	    $type=0;
	}
	print FH "$n $x $y $type\n";
	$n++;
    }
}


$x_del=6.22;
$y_del=4-31.9+2.5;
for ($j=0;$j<16;$j++) {
    for ($i=0;$i<16;$i++) {
	$x=$x_del+0.5*$j;
	$y=$y_del+0.5*$i;
	$type=1;
	print FH "$n $x $y $type\n";
	$n++;
    }
}


close(FH);

exit;
