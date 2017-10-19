#!/usr/bin/perl


open(DIR,"ls *.rss.fits |");
while($file=<DIR>) {
    chop($file);
    if (($file !~ "mask")&&($file !~ "mos")&&($file !~ "rad")&&($file !~ "std")) {
	print "$file\n";
	$cube=$file;
	$cube =~ s/rss/cube/;
	$B=$file;
	$B =~ s/rss/B/;
	$call="cp ".$file." mos.rss.fits";
	system($call);

	$call="rss2rss_mask_ppak.pl mos.rss.fits mos.pt.txt mask.rss.fits ".$file." 1";
	system($call);

	$call="rm -r ".$cube." ";
	system($call);
	
	$call="rss2cube.tcl ".$file." mos.pt.txt ".$cube." 3 1e-12 1";
	system($call);
	$call="get_slice.pl ".$cube." img /disk-b/sanchez/ppak/legacy/slice_B.conf";
	system($call);
	$call="cp img_B_3900_4550.fits ".$B." ";
	system($call);


    }

}
close(DIR);
