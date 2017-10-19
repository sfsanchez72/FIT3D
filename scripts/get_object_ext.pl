#!/usr/bin/perl


open(DIR,"ls run*b.fits |");
print "# file object date ra dec AZ EL AIRMASS JULIAN_DATE\n";

while ($file=<DIR>) {
    chop($file);

    $call="read_img_header.pl ".$file." OBJECT > get_object_ext.out";
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	$out =~ s/ //g;
	$out =~ s/$file//g;
	$name=$file;
	$object=$out;
#	($name,$object)=split(" ",$out);
#	$object =~ s/ //g;

    }
    close(CALL);

    $call="read_img_header.pl ".$file." RA > get_object_ext.out";
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	($name,$ra)=split(" ",$out);
    }
    close(CALL);

    $call="read_img_header.pl ".$file." DEC > get_object_ext.out";
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	($name,$dec)=split(" ",$out);
    }
    close(CALL);

    $call="read_img_header.pl ".$file." DATE > get_object_ext.out";
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	($name,$date)=split(" ",$out);
    }
    close(CALL);


    $call="read_img_header.pl ".$file." AIRMASS > get_object_ext.out";
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	($name,$airmass)=split(" ",$out);
    }
    close(CALL);

    $call="read_img_header.pl ".$file." JUL-DATE > get_object_ext.out";
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	($name,$jul)=split(" ",$out);
    }
    close(CALL);


    $call="read_img_header.pl ".$file." HIERARCH | grep AZ_START > get_object_ext.out"; 
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	@data=split(" ",$out);
	$AZ=$data[5];
    }
    close(CALL);

    $call="read_img_header.pl ".$file." HIERARCH | grep EL_START > get_object_ext.out"; 
    system($call);
    open(CALL,"<get_object_ext.out");
    while($out=<CALL>) {
	chop($out);
	@data=split(" ",$out);
	$EL=$data[5];
    }
    close(CALL);




    print "$file $object $date $ra $dec $AZ $EL $airmass $jul\n";

}
exit;
