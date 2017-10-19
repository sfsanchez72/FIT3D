#!/usr/bin/perl
#
#
#

if ($#ARGV<4) {
    print "USE:velocity.pl FITSFILE LINE_FILE START_WAVELENGTH END_WAVELENGTH VIS\n";    
    exit;
}

$fitsfile=$ARGV[0];
$em_file=$ARGV[1];
$w_start=$ARGV[2];
$w_end=$ARGV[3];
$vis=$ARGV[4];

if (($vis!=0)&&($vis!=1)) {
    $vis=1;
}

$n=0;
open(FH,"<$em_file");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $w[$n]=$data[0];
    $name[$n]=$data[1];
    $n++;
}
close(FH);
print "$n lines found\n\n";



print "Number of systems to include:";
$n_s=<stdin>;
if ($n_s!=int($n_s)) {
    print "\nError\n";
    exit;
}

for ($i=0;$i<$n_s;$i++) {
    print "Redshift of the first system:";
    $z[$i]=<stdin>;
    chop($z[$i]);
    print "Flexibility on the Redshift (%):";
    $flex_w=<stdin>;
    chop($flex_w);
    print "Velocity dispersion (Amstrongs!) (vel,min,max):";
    $tmp=<stdin>;
    chop($tmp);
    ($dz[$i],$mdz[$i],$Mdz[$i])=split(" ",$tmp);        
    print "Flux level of the first system:";
    $f[$i]=<stdin>;
    chop($f[$i]);
    print "MIN Flux level of the first system:";
    $f_m[$i]=<stdin>;
    chop($f_m[$i]);
    print "MAX Flux level of the first system:";
    $f_M[$i]=<stdin>;
    chop($f_M[$i]);

}

$n_line;
for ($i=0;$i<$n_s;$i++) {
    $nl[$i]=0;
    for ($j=0;$j<$n;$j++) {
	$w1=$w[$j]*(1+$z[$i]);
	if (($w1>$w_start)&&($w1<$w_end)) {
	    print "$name[$j] $w[$j] at $w1 found\n";
	    print "Include? (y/n):";
	    $include=<stdin>;
	    chop($include);
	    if ($include !~ "n") {
		$wl[$nl[$i]][$i]=$w1;
		print "Flux (Flux, min, Max)? (y/n):";
		$tmp=<stdin>;
		chop($tmp);
		($flux[$nl[$i]][$i],$min_flux[$nl[$i]][$i],$max_flux[$nl[$i]][$i])=split(" ",$tmp);
		$nl[$i]++;
		$n_line++;	
		print "Included\n";
	    }	
	}    
    }
    print "System $i contains $nl[$i] lines\n";
}
print "The total number of lines is $n_line\n";
$n_mod=$n_line+1;




$flex_dw=0.05;

open(FH,">tmp.config");
print FH "0 $n_mod 0.2 0.001\n";
for ($i=0;$i<$n_s;$i++) {
    $j=0;
    $a=$wl[$j][$i];
    $p_a=$wl[$j][$i]-$flex_w*$wl[$j][$i]-$dz[$i];
    $p_a2=$wl[$j][$i]+$flex_w*$wl[$j][$i]+$dz[$i];
    $d_a=$dz[$i];
    $p_da=$mdz[$i];
    $p_da2=$Mdz[$i];
    $ff=$flux[$j][$i];
    $mff=$min_flux[$j][$i];
    $Mff=$max_flux[$j][$i];
    print FH "gauss1d\n";
    print FH "$a\t 1\t $p_a\t $p_a2\t -1\n"; 
    print FH "$ff\t 1\t $mff\t $Mff\t -1\n"; 
    print FH "$d_a\t 1\t $p_da\t $p_da2\t -1\n"; 
    for ($k=0;$k<6;$k++) {
	printf FH "0\t 0\t 0\t 0\t -1\n";
    }
    for ($j=1;$j<$nl[$i];$j++) {
	$a=$wl[0][$i];
	$p_a=$wl[$j][$i]-$wl[0][$i];
	$p_a2=0.0;
	$d_a=$dz[$i];
	$p_da=$mdz[$i];
	$p_da2=$Mdz[$i];
	$ff=$flux[$j][$i];
	$mff=$min_flux[$j][$i];
	$Mff=$max_flux[$j][$i];
	print FH "gauss1d\n";
	print FH "$a\t 1\t $p_a\t $p_a2\t 1\n"; 
	print FH "$ff\t 1\t $mff\t $Mff\t -1\n"; 
	print FH "$d_a\t 1\t $p_da\t $p_da2\t 1\n"; 
	for ($k=0;$k<6;$k++) {
	    printf FH "0\t 0\t 0\t 0\t -1\n";
	}
    }
}

print "BACKGROUND order:";
$ord_back=<stdin>;
chop($ord_back);
if ($ord_back>8) { 
    $ord_back=8;
}
if ($ord_back<1) {
    $ord_back=1;
}

print "BACKGROUND LEVEL:";
$back=<stdin>;
chop($back);
$back_pa=0.1*$back;
$back_pa2=10**$back;

print FH "poly1d\n";
print FH "$back\t 1\t $back_pa\t $back_pa2\t -1\n"; 
for ($k=0;$k<$ord_back;$k++) {
    printf FH "0.0001\t 1\t 0.0\t 5.0\t -1\n";
}
for ($k=$ord_back;$k<(8-$ord_back-1);$k++) {
    printf FH "0\t 0\t 0\t 0\t -1\n";
}

close(FH);

system("rm out.model");
system("rm model.fits");
system("rm res.fits");
$call="velocity -if ".$fitsfile." -cf tmp.config -vis ".$vis." -flex 0.2";
system($call);


#
# Once done the analysis, we need to create
# the position Intensity, Velocity and Dispersion Maps.
#




exit;


