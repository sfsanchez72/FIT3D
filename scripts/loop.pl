#!/usr/bin/perl


$file="/disk-b/sanchez/ppak/legacy/templates/list_of_list.txt";

$nf=0;

open(FILE,"<$file");
while($file=<FILE>) {
    chop($file);
    if ($file !~ "#") {
	$files[$nf]=$file;
	$nn[$nf]=0;
	$k=0;
	open(NAME,"<$file");
	while($name=<NAME>) {
	    chop($name);
	    if ($name !~ "#") {
		$nam_file[$nf][$k]=$name;
		$k++;
	    }
	}
	close(NAME);
#	print "$nf $k\n";
	$nn[$nf]=$k;
       	$nf++;
    }
}
close(FILE);

print "nf=$nf\n";
$NC=0;
$name_now=loop("",0);
#print "name_now = $name_now";


for ($j=0;$j<$NC;$j++) {
    print "$j/$NC $LIST[$j]\n";
}

exit;


sub loop {
    my $name=$_[0];
    my $NF=$_[1];
    my $i=0;
    if ($NF<$nf) {
	for ($i=0;$i<$nn[$NF];$i++) {
	    my $name_now=$name." ".$nam_file[$NF][$i];
	    my $val=loop($name_now,$NF+1);
	}
    } else {
	$LIST[$NC]=$name;
	$NC++;
#	print "$name\n";
	return $name;
    }
}

