#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");

if ($#ARGV<4) {
    print "USE: trigger_cal.pl INPUT OBJECT CAL SKY PLOT\n";
    exit;
}

$inputfile=$ARGV[0];
$objfile=$ARGV[1];
$calfile=$ARGV[2];
$skyfile=$ARGV[3];
$spy=$ARGV[4];

$n=0;
$nc=15;
$ns=36;
open(FH,"<ppak_positions.txt");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $type[$n]=0;
    if (($data[1]>400)&&($data[1]<500))    {
	$type[$n]=1;
#	print "\n";
    } 

    if (($data[1]>500)&&($data[1]<900))    {
	$type[$n]=2;
#	print "\n";
    }

    if ($data[1]<400)    {
#	print "$data[0] ";
    }


    if ($data[1]<900) {
	$n++;
    }

    
}
close(FH);

#print "N=$n\n";

print "Reading input file\n";
$naxes=read_naxes($inputfile);
($n1,$n2) = @$naxes;
@input=read_img($inputfile);
print "Done\n";
print "$n1 $n2\n";


$jc=0;
$js=0;
$ji=0;
print "triggering\n";
for ($j=0;$j<$n2;$j++) {
$y_max=-100000;
$y_min=100000;

#    print "type=$type[$j] $jc/15 $js/36\n";
    if ($type[$j]==2) {
	for ($i=0;$i<$n1;$i++) {
	    $cal[$jc][$i]=$input[$j][$i];
	    $y[$i]=$input[$j][$i];
	    $x[$i]=$i;
	    if ($y_max<$y[$i]) {
		$y_max=$y[$i];
	    }
	    if ($y_min>$y[$i]) {
		$y_min=$y[$i];
	    }
	}
	if ($spy==1) {
	    pgbegin(0,"/xs",1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.6);           # Set character height
	    pgenv(0,$n1,$y_min,$y_max,0,0);
	    pgsch(1.6);           # Set character height
	    pglabel("Id.","Counts","Cal $jc");
	    pgsch(2.2);           # Set character height
	    pgsci(1);
	    pgline($n1,\@x,\@y);    
	    pgsci(1);
	    pgclose;
	    pgend;
	    <stdin>;
	}


	$jc++;

    }

    if ($type[$j]==1) {
	for ($i=0;$i<$n1;$i++) {
	    $sky[$js][$i]=$input[$j][$i];
	    $y[$i]=$input[$j][$i];
	    $x[$i]=$i;
	    if ($y_max<$y[$i]) {
		$y_max=$y[$i];
	    }
	    if ($y_min>$y[$i]) {
		$y_min=$y[$i];
	    }
	}


	if ($spy==1) {
	    pgbegin(0,"/xs",1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.6);           # Set character height
	    pgenv(0,$n1,$y_min,$y_max,0,0);
	    pgsch(1.6);           # Set character height
	    pglabel("Id.","Counts","Sky $js");
	    pgsch(2.2);           # Set character height
	    pgsci(1);
	    pgline($n1,\@x,\@y);    
	    pgsci(1);
	    pgclose;
	    pgend;
	    <stdin>;
	}


	    $js++;

    }

    if ($type[$j]==0) {
	for ($i=0;$i<$n1;$i++) {
	    $obj[$ji][$i]=$input[$j][$i];
	    $y[$i]=$input[$j][$i];
	    $x[$i]=$i;
	    if ($y_max<$y[$i]) {
		$y_max=$y[$i];
	    }
	    if ($y_min>$y[$i]) {
		$y_min=$y[$i];
	    }
	}
	if ($spy==1) {
	    pgbegin(0,"/xs",1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.6);           # Set character height
	    pgenv(0,$n1,$y_min,$y_max,0,0);
	    pgsch(1.6);           # Set character height
	    pglabel("Id.","Counts","IFU $ji");
	    pgsch(2.2);           # Set character height
	    pgsci(1);
	    pgline($n1,\@x,\@y);    
	    pgsci(1);
	    pgclose;
	    pgend;
	    <stdin>;
	}


	$ji++;

    }


#    print "$j $jc $js $type[$j]\n";
}

print "Writting the results\n";
system("rm $calfile");
write_img($calfile,$n1,15,\@cal);
print "$calfile saved\n";
system("rm $skyfile");
write_img($skyfile,$n1,36,\@sky);
print "$skyfile saved\n";
system("rm $objfile");
write_img($objfile,$n1,331,\@obj);
print "$objfile saved\n";

exit;
