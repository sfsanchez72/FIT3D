#!/usr/bin/perl


if ($#ARGV<3) {
    print "change_line.pl input_line output_line inputfile outputfile\n";
    exit;
}

$inline=$ARGV[0];
$outline=$ARGV[1];
$infile=$ARGV[2];
$outfile=$ARGV[3];
open(FH,"<$infile");
open(FHOUT,">$outfile");
print "$inline\n";
print "$outline\n";
while($line=<FH>) {
    chop($line);
    $line =~ s/$inline/$outline/g;
    print FHOUT "$line\n";
}
close(FHOUT);
close(FH);

exit;
