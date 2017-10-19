#!/usr/bin/perl
#
use Statistics::OLS;
use Math::FFT;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);
use Math::Approx;
use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<0) {
    $dev="";
} else {
    $dev=$ARGV[0];
}
$infile="object.txt";
$n=0;
open(FH,"<$infile");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$file[$n]=$data[0];
	$obj[$n]=$data[1];
	$fnow=$data[0];
	$fnow =~ s/b//;
	open(HEADER,"dump_header.pl $fnow |");
	while($header=<HEADER>) {
	    chop($header);
	    @DATA=split(/\|/,$header);
	    if (($DATA[0] =~ "EXPTIME")&&($DATA[0] !~ "EXPTIMEM")) {
		$exptime{$fnow}=$DATA[1];
	    }
	    if ($DATA[0] =~ "DATE_OBS") {
		$date{$fnow}=$DATA[1];
	    }
	    if ($DATA[0] =~ "OBSERVER") {
		$observer=$DATA[1];
	    }
	}
	close(HEADER);
	$DATE=$date{$fnow};
        ($day,$ut)=split("T",$DATE);
	($h,$m,$s)=split(":",$ut);
	if ($n==0) {
	    $DAY=$day;
	    $offset=-24;
	} else {
	    if ($day ne $day_last) {
		$offset=0;
	    }
	}
	$h=$h+$offset;
	$t=$h+$m/60+$s/3600;
	$time{$fnow}=$t;
	$day_last=$day;

	$n++;
    }
}
close(FH);

$DAY =~ s/ //g;
if ($dev eq "") {
    $dev="cal_obs_seq.".$DAY.".ps/CPS";
}

$Tu=0;
$X_end=0;
$X_start=0;
pgbegin(0,$dev,1,1);
pgsfs(1.2);
pgscf(2);             # Set character font
pgslw(2);             # Set line width
pgsch(1.4);           # Set character height
pgenv(-9,9,-1,3,0,0);
pglabel("U.T.","Observation","$DAY");
for ($i=0;$i<$n;$i++) {
    $fnow=$file[$i];
    $obj_now=$obj[$i];
    $X1=$time{$fnow};
    $X2=$time{$fnow}+$exptime{$fnow}/3600;
    if ($obj_now =~ "obj") {
	$Y=2;
	$color=1;
	$X_end2=$X2;
	if (($obj_now =~ "cont")||($obj_now =~ "arc")) {
	    $Y=1.5;
	    $color=4;
	} else {
	    $Tu=$Tu+$exptime{$fnow}/3600;
	}
    } else {
	$Y=0;
	$color=2;
	if ($obj_now =~ "sky") {
	    $Y=1;
	    $color=8;
	    if ($obj_now =~ "sky_1") {
		$X_start=$X2;
	    } 
	    if (($obj_now =~ "sky_2")&&($X_end==0)) {
		$X_end=$X1;
	    }
	}
    }
    $Y1=$Y-0.25;
    $Y2=$Y+0.25;
#    print "$fnow $obj_now $color\n";
    
    pgsci($color);
    pgrect($X1,$X2,$Y1,$Y2);
}

if ($X_end==0) {
    $X_end=$X_end2;
}
$nl=$X_end-$X_start;
$observer =~ s/ //g;
$nl=apr($nl);
#print "# (1) Observer\n";
#print "# (2) Night Length (hours)\n";
#print "# (3) Hours on target (hours)\n";
print "$observer,$nl,$Tu\n";
pgsci(1);
pgsch(1.2);
$text=$observer.", Night Length: ".$nl." h., On Target:".$Tu." h.";
pgptxt(-8.5,2.7,0,0,"$text");


pgsci(1);
pgclose;
pgend;


exit;





