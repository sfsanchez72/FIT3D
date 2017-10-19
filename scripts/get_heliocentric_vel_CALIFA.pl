#!/usr/bin/perl

if ($#ARGV<0) {
    print "USE: get_heliocentric_vel_PPAK.pl FITSFILE\n";
    exit;
}

$file=$ARGV[0];

open(FH,"dump_header.pl $file PPAK |");
while($line=<FH>) {
    chop($line);
    ($key,$val)=split(/\|/,$line);
    $key =~ s/ //g;
    $val =~ s/ //g;
    if ($key eq "PPAK_PPAK_P1_F0_RA") {
	$RA=$val;
    }
    if ($key eq "PPAK_PPAK_P1_F0_DEC") {
	$DEC=$val;
    }
    if ($key eq "PPAK_PPAK_P1_F0_DATE_OBS") {
	$DATE=$val;
    }

}
close(FH);
#print "$RA |  $DEC | $DATE"; exit;

$rh=($RA/15);
$RH=int($rh);
if ($RH<10) {
    $RH="0".$RH;
}
$rm=($rh-int($rh))*60;
$RM=int($rm);
if ($RM<10) {
    $RM="0".$RM;
}
$rs=($rm-int($rm))*60;
$RS=int($rs);
if ($RS<10) {
    $RS="0".$RS;
}
$ra=$RH.":".$RM.":".$RS;

if ($DEC>0) {
    $sig="+";
} else {
    $sig="-";
}
$dh=abs($DEC);

$DH=int($dh);
if ($DH<10) {
    $DH="0".$DH;
}
$dm=($dh-int($dh))*60;
$DM=int($dm);
if ($DM<10) {
    $DM="0".$DM;
}
$ds=($dm-int($dm))*60;
$DS=int($ds);
if ($DS<10) {
    $DS="0".$DS;
}
$dec=$sig.$DH.":".$DM.":".$DS;

@dat=split(/\-/,$DATE);
$eq="J".$dat[0];
$eq =~ s/ //g;

#$eq="J2000";
#print "$ra $dec $DATE $eq\n";
$call="get_heliocentric_vel.pl CAHA ".$ra." ".$dec." ".$eq." ".$DATE;
system($call);



exit;

