#!/usr/bin/perl


if ($#ARGV<1) {
    print "change_line.pl input_line output_line\n";
    exit;
}

$inline=$ARGV[0];
$outline=$ARGV[1];

open(DIR,"ls * |");
while($file=<DIR>) {
    chop($file);
    if ($file ne "change_line.pl") {
	system("cp $file tmp");
	open(FH,"<tmp");
	open(FHOUT,">$file");
	while($line=<FH>) {
	    chop($line);
	    $line =~ s/$inline/$outline/;
	    print FHOUT "$line\n";
	}
	close(FHOUT);
	close(FH);
	system("rm tmp");
    }
}
close(DIR);
exit;
