#!/usr/bin/perl


if ($#ARGV<4) {
    print "USE: reduce_all.pl INFILE TRACING_REF WAVELENGTH_REF EXPTIME CONFIGFILE AIRMASS EXTINCTION SIGMA\n";
    exit;
}

$infile=$ARGV[0];
$trace=$ARGV[1];
$wave=$ARGV[2];
$exptime=$ARGV[3];
#$cdelt=$ARGV[4];
$config=$ARGV[4];
$airmass=$ARGV[5];
$extinction=$ARGV[6];
$sigma=$ARGV[7];
$infile =~ s/.fits//;
$trace =~ s/.fits//;
$wave =~ s/.fits//;

$prefix=$infile;
$prefix =~ s/.fits//;
$log="red_".$prefix.".log";



$t_start=time;
open(FH,"<$config");
open(OUT,">$log");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	$line =~ s/PREFIX/$prefix/g;
	$line =~ s/TRACE/$trace/g;
	$line =~ s/CALIBW/$wave/g;
	$line =~ s/INFILE/$infile/g;
	$line =~ s/EXPTIME/$exptime/g;
	$line =~ s/AIRMASS/$airmass/g;
	$line =~ s/EXTINCTION/$extinction/g;
	$line =~ s/SIGMA/$sigma/g;
#	$line =~ s/CDELT/$cdelt/g;
	print "$line\n";
	print OUT "$line\n";
	system($line);
	if ($line =~ "peak_find_width") {
	    open(W,"<tjunk.mean.width");
	    $w=<W>;
	    @dat=split(" ",$w);
	    $sigma=$dat[3];
	    close(W);
	}

	print "DONE\n";
    }
}


close(OUT);
close(FH);
$delta_t=(time-$t_start)/60;
print "Consumed time: $delta_t\n";
exit;


