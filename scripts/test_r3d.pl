#!/usr/bin/perl

open(DIR,"ls *.pl|");
while($file=<DIR>) {
    chop($file);
    if (($file ne "my.pl")&&($file ne "test_r3d.pl")) {
	$ex="./".$file;
	system($ex);
    }
}
close(DIR);

exit;
