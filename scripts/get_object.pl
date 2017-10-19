#!/usr/bin/perl


open(DIR,"ls run*b.fits |");
while($file=<DIR>) {
    chop($file);
    $call="read_img_header.pl ".$file." OBJECT";
    system($call);
}
close(DIR);

exit;
