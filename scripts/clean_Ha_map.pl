#!/usr/bin/perl
#!/c/strawberry/perl/bin/perl.exe
#
use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;

#require "/work1/ssa/perl/MY/my.pl";
#require 'c:\mingw\msys\1.0\home\sanchez\sda2\code\R3D\my.pl';
require "/home/sanchez/sda2/code/R3D/my.pl";


if ($#ARGV<1) {
    print "USE: imarith.pl INPUT1.FITS INPUT2.FITS [MAX_VEL] [MIN_VEL]\n";
    exit;
}

$infile1=$ARGV[0];
$infile2=$ARGV[1];
$max_vel=1e12;
if ($#ARGV==2) {
  $max_vel=$ARGV[2];
}

$min_vel=0;
if ($#ARGV==3) {
  $max_vel=$ARGV[2];
  $min_vel=$ARGV[3];
}


$a_in1=rfits($infile1);
$h=$a_in1->gethdr;
$a_in2=rfits($infile2);
($nx,$ny)=$a_in1->dims();
for ($i=0;$i<$nx;$i++) {
  for ($j=0;$j<$ny;$j++) {
    $vel=$a_in1->at($i,$j);
    $mask=$a_in2->at($i,$j);
    if ($mask==1) {
      $a[$k]=$vel;
      $k++;
    }
  }
}

$pdl_a=pdl(@a);
@stats=stats($pdl_a);
#print "@stats\n"; exit;
for ($i=0;$i<$nx;$i++) {
  for ($j=0;$j<$ny;$j++) {
    $vel=$a_in1->at($i,$j);
    $mask=$a_in2->at($i,$j);
    if ((($vel-$stats[2])>300)||($vel>$max_vel)) {
      $lo=6562*(1+$vel/300000);
      $new_vel=(($lo/6583)-1)*300000;
      set($a_in1,$i,$j,$new_vel);
      set($a_in2,$i,$j,1);
    }
    if ($vel<$min_vel) {
      set($a_in2,$i,$j,0);
    }
  }
}


$a_in1->wfits($infile1);
$a_in2->wfits($infile2);

exit;
