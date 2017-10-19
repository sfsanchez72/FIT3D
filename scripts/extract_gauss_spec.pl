#!/usr/bin/perl
#
#
#
use PGPLOT;


use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");

# This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s the and output file
# like the one created by "fit_3D_fg.pl"
# and fix one parameter to do
# a force fitting latter
if ($#ARGV<4) {
    print "USE: extract_moffat_spec.pl MODEL INPUT_FILE OUT_FILE CRVAL CDELT\n";        
    exit;
}

#$n_param=$ARGV[0];
$n_model=$ARGV[0];
$input_file=$ARGV[1];
$out_file=$ARGV[2];
$crval=$ARGV[3];
$cdelt=$ARGV[4];

#
# Reading the input file
#
$y_min=1000000;
$y_max=-1000000;
$n=0;
open(INPUT,"<$input_file");
while($input=<INPUT>) {
    chop($input);
    $n_mod=$input;
    for ($i=0;$i<$n_mod;$i++) {
	$input=<INPUT>; chop($input);
	@data=split(" ",$input);
	$model[$n]=$data[0];
	for ($j=0;$j<9;$j++) {
	    $a[$i][$j][$n]=$data[1+2*$j];
	    $ia[$i][$j][$n]=$data[1+2*$j+1];
	}
#
#
    if ($i==$n_model) {
	$sigma=$a[$i][2][$n];
	$rd=$a[$i][3][$n];
	$a_fix[$n]=2*3.1416*$sigma*($rd**2);
	$s_y_out2[$n]=$a_fix[$n];
    }



    }
    $x[$n]=$n;
    $n++;
};
close(INPUT);



open(OUTPUT,">$out_file");
for ($k=0;$k<$n;$k++) {
    $wave=$crval+$cdelt*$k;
    print OUTPUT "$k $wave $a_fix[$k]\n";
};
close(OUTPUT);

exit;
#




sub mean {
  local(@data)=@_;
  my $sum=0; 
  my $j;
  for ($j=0;$j<$#data;$j++) {
      $sum=$sum+$data[$j];
  }
  my $mean;
  if ($#data>0) {
     $mean = $sum/($#data);
  } else {
     $mean=-666;
  }
  return $mean;
}


sub stats {
  local(@data)=@_;
  my $mean = mean(@data);
  my $stddev = 0;
  my $j;
  my $sum=0; 
  for ($j=0;$j<$#data;$j++) {
      $sum=$sum+($data[$j]-$mean)**2;
  }
  if ($#data>0) {
      if ($#data>1) {
	  $sum=$sum/($#data-1);
      } else {
	  $sum=$sum/($#data);
      } 
      $stddev=sqrt($sum);
  } else {
      $stddev=-666;
  }
  my @out=($mean,$stddev);

  return @out;
}



