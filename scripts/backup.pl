#!/usr/bin/perl


@DIR=("/home/sanchez/sda2/","/home/sanchez/sda1/bin/","/home/sanchez/sda1/perl/","/home/sanchez/Documents/","/home/sanchez/Mail/","/disk-a/sanchez/","/disk-b/sanchez/","/disk-c/sanchez/");

foreach $dir (@DIR) { 
    @SEC=split(/\//,$dir);
    my $new="/media/usb1/sanchez";
    foreach $sec (@SEC) {
	$new=$new."/".$sec;
	$call="mkdir ".$new;
	system($call);
    }

#    $out="/media/disk/sanchez".$dir;
#    $out="/media/35193363-50c1-4bfe-9c56-7b798354dca5/sanchez".$dir;
    $out="/media/usb1/sanchez".$dir;
    $call="cp -r -u ".$dir."* ".$out;
    system($call);    
    print "'$dir' in the backup disk\n";
}

exit;




