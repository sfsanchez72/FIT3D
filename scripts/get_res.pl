#!/usr/bin/perl


if ($#ARGV<0) {
    print "get_res.pl GRATING_NAME\n";
    exit;
}
$GRAT=$ARGV[0];
$CONF="/disk-b/sanchez/ppak/legacy/";

$call="rm -rf get_res.out";
system($call);
print "We use the spectra: spec_*.msky.txt as reference\n";
open(DIR,"ls spec_*.msky.txt | ");
while($line=<DIR>) {    
    chop($line);
    @data=split(/\_/,$line);
    $name=$data[1];
    if (($line !~ "obj")&&($name ne "")) {
#	print "$line\n";
	$call="fit_spec_back.pl ".$line." ".$CONF."/emission_lines_fit.sky 3700 7200 none 0 ".$CONF."/mask_emission_lines_fit.sky ".$CONF."/SKY_".$GRAT.".config /null > junk.junk.junk";
	system($call);
	$call="cat out.fit_spectra | grep eline >> get_res.out";
	system($call);
    }

}
close(DIR);

$call="table_hist.pl get_res.out 5 get_res.ps/CPS 10 > med_get_res.out";
system($call);

exit;

