#!/usr/bin/perl


open(DIR,"ls run*_?????.fits |");
while($file=<DIR>) {
    chop($file);
    $call="read_img_header.pl ".$file." OBJECT > line.txt";
    system($call);
    open(LINE,"<line.txt");
    $line=<LINE>;
    close(LINE);
    chop($line);
    $line =~ s/ //g;
    $line =~ s/$file//g;
    print "$file  $line\n";
}
close(DIR);

exit;
