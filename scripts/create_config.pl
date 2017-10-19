#!/usr/bin/perl

if ($#ARGV<1) {

    print "USE: create_config.pl ARC.disp.id PREFIX\n";
    exit;
}


$input=$ARGV[0];
$prefix=$ARGV[1];

open(FH,"<$input");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$w=int($data[1]);
	$output=$prefix.".".$w.".config";
	open(OUT,">$output");
	print OUT "0 2 0.2 0.001\n";
	print OUT "eline\n";
	print OUT "$data[1]	 0	 0	 0	 -1\n";
	print OUT "1e6	 1	 0	 1e12	 -1\n";
	print OUT "1.5	 1	 0.1	 30	 -1\n";
	print OUT "0	 1	 -500	 500	 -1\n";
	print OUT "0	 0	 0	 0	 -1\n";
	print OUT "0	 0	 0	 0	 -1\n";
	print OUT "0	 0	 0	 0	 -1\n";
	print OUT "0	 0	 0	 0	 -1\n";
	print OUT "0	 0	 0	 0	 -1\n";
	print OUT "poly1d\n";
	print OUT "1       1       -1e13   1e13    -1\n";
	print OUT "0       0       -1e13   1e13    -1\n";
	print OUT "0        0       -1e13   1e13    -1\n";
	print OUT "0        0       -1e13   1e13    -1\n";
	print OUT "0        0       -1e13   1e13    -1\n";
	print OUT "0        0       -1e13   1e13    -1\n";
	print OUT "0        0       0       0       -1\n";
	print OUT "0        0       0       0       -1\n";
	print OUT "0        0       0       0       -1\n";
	close(OUT);

    }    
}
close(FH);
