#!/usr/bin/perl


if ($#ARGV<1) {
    print "USE: spots_list.pl objects.txt spots_list [DX DY]\n";
    exit;
}

$infile=$ARGV[0];
$spots_file=$ARGV[1];
$DX=$ARGV[2];
$DY=$ARGV[3];

system("rm -f *.spots");

open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$name=$data[0];
	$type=$data[1];
	if (((($type =~ "arc")||($type =~ "obj")))&&($type !~ "cont")) {
	    $name_out=$name;
	    $name_out =~ s/fits/spots/;
	    $call="spots_center.pl ".$name." ".$spots_file." 3 3 ".$name_out." ".$DX." ".$DY;
	    print "$call\n";
	    system($call);
	    print "$name processed\n";
	}

    }
}

exit;
