#!/usr/bin/perl
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

if ($#ARGV<3) {
    print "get_resolution_map.pl ARC.disp_cor.fits ID.file output.file DEV\n";
    exit;
}

$arc_file=$ARGV[0];
$id_file=$ARGV[1];
$out_file=$ARGV[2];
$dev=$ARGV[3];

$na=0;
open(FH,"<$id_file");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	@data=split(" ",$line);
	$wa[$na]=$data[1];
	$na++;
    }
    
}
close(FH);

#print "NA=$na\n";
$call="rm -f ".$out_file;
system($call);
$call="rm -f get_resolution_map.log";
system($call);
$pdl=rfits($arc_file);
($nx,$ny)=$pdl->dims;
$crval=$pdl->hdr->{"CRVAL1"};
$cdelt=$pdl->hdr->{"CDELT1"};
$crpix=$pdl->hdr->{"CRPIX1"};
if ($cdelt==0) {
    $cdelt=1;
}

for ($j=0;$j<$ny;$j++)  {
    open(FH,">spec_tmp.txt");
    for ($i=0;$i<$nx;$i++) {
	$k=$i+1;
	$wave[$i]=$crval+$cdelt*($i+1-$crpix);
	$flux[$i]=$pdl->at($i,$j);
	if ($flux[$i] eq "BAD") {
	    $flux[$i]=0;
	}
	print FH "$k $wave[$i] $flux[$i]\n";
    }
    close(FH);


#    $call="img2spec.pl ".$arc_file." ".$j." spec_tmp.txt";
#    system($call);
    $call="echo '$na' >> ".$out_file;
    system($call);
    for ($i=0;$i<$na;$i++) {
      open(TMP,">tmp.arc.config");
      print TMP "0 2 1 0.5\n";
      print TMP "eline\n";
      print TMP "$wa[$i]	 0	 0	 0	 -1\n";
      print TMP "1e4	 1	 0	 1e20	 -1\n";
      print TMP "2.5	 1	 0.1	 20	 -1\n";
      print TMP "0	 1	 -100	 100	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
      print TMP "poly1d\n";
      print TMP "0.1	 1	 -1e13	 1e13	 -1\n";
      print TMP "0.01	 1	 -1e13	 1e13	 -1\n";
      print TMP "0	 0 	 -1e13	 1e13	 -1\n";
      print TMP "0	 0	 -1e13	 1e13	 -1\n";
      print TMP "0	 0	 -1e13	 1e13	 -1\n";
      print TMP "0	 0	 -1e13	 1e13	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
      print TMP "0	 0	 0	 0	 -1\n";
   close(TMP);
      $wa_min=$wa[$i]-25;
      $wa_max=$wa[$i]+25;
      $call="fit_spec_back.pl spec_tmp.txt none ".$wa_min." ".$wa_max." none 0 none tmp.arc.config ".$dev." >> get_resolution_map.log";
      system($call);
      $call="cat out.fit_spectra | grep eline >> ".$out_file;
      system($call);
      print "$i/$na $j/$ny\n";
#      print "Press Enter\n"; <stdin>;
    }
}


  exit;

#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
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


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";

if ($#ARGV<1) {
    print "USE: img2spec.pl INPUT_FILE.FITS NY OUTPUTFILE.txt\n";
    exit;
}

$input=$ARGV[0];
$NY=$ARGV[1];
$output=$ARGV[2];


$y_min=1e12;
$y_max=-1e12;
$n=0;
$pdl=rfits("$input");
($nx,$ny)=$pdl->dims;

exit;
