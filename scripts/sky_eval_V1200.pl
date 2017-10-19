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
$DIR_CONF="/disk-b/sanchez/ppak/legacy/";

if ($#ARGV<0) {
    print "USE: auto_redu.pl DATE(YYYYMMDD) GLUE[0=NO/1=YES] [DX DY spots]\n";
    exit;
}

#$call="get_object_ext.pl > object_ext.txt";
#system($call);

$nno=0;
open(FH,"<object_ext.txt");
while($line=<FH>) {
    chop($line);
    $a_line[$nno]=$line;
    $nno++;
}
close(FH);


#
# SKYs
#
open(SKY,"ls -ltra spec_*msky.txt |");
open(OUTSKY,">auto_redu_V1200.sky_list");
open(OUTMAG,">auto_redu_V1200.sky_mag");
print OUTMAG "# sky_file V_mag_sky\n";
while($line=<SKY>) {
    chop($line);
    @data=split(" ",$line);
    $file=$data[7];
    $nam=$file;
    $cut="spec_";
    $nam =~ s/$cut//;
    $cut=".msky.txt";
    $nam =~ s/$cut//;
    ($pre,$point)=split(/\_/,$nam);


    $NAME=$point."_".$pre;
    if (($file !~ "run")&&($file !~ "obj")&&($file !~ "vstar")) {
	print OUTSKY "$file\n";
	$call="mag_filter.pl ".$DIR_CONF."/filters/B_short.txt 19.47 ".$file." 0 > mag_filter.out";
	system($call);
	print "$call\n";
	open(MAG,"<mag_filter.out");
	$a_mag=<MAG>;	
	chop($a_mag);
	$line_now="";
	for ($i=0;$i<$nno;$i++) {
	    @DATA=split(" ",$a_line[$i]);
	    $NAM=$DATA[1];
	    $NAME =~ s/\+//;
	    $NAM =~ s/\+//;
	    if ($NAM =~ $NAME) {
		print "$NAM | $NAME\n";

		$line_now=$a_line[$i];
	    }
	}

	print OUTMAG "$a_mag $line_now\n";

	$call="get_flux_stats_cube.pl ".$pre.".cube.fits 4500 40 3400 0.7 SN_".$pre.".ps/CPS 1 > SN_".$pre.".out";
	system($call);
    }
}
close(OUTMAG);
close(OUTSKY);
close(SKY);


$call="spectra_plot.pl auto_redu_V1200.sky_list skys_plot.ps/CPS -1 10 3700 4800";
system($call);

$label1="Julian DATE";
$label2="Night-sky B-band (mag/arcsec\u2\d)";
$call="table_plot.pl auto_redu_V1200.sky_mag 11 1 '".$label1."' '".$label2."' sky_mag_date.ps/CPS";
system($call);
print "$call";


$label1="Airmass";
$label2="Night-sky B-band (mag/arcsec\u2\d)";
$call="table_plot.pl auto_redu_V1200.sky_mag 10 1 '".$label1."' '".$label2."' sky_mag_airmass.ps/CPS";
system($call);

exit;
