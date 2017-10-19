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
if ($#ARGV<5) {
    print "USE: mean_force.pl PARAM MODEL N_SIGMA BOX_FWHM INPUT_FILE OUT_FILE\n";        
    exit;
}

$n_param=$ARGV[0];
$n_model=$ARGV[1];
$n_sigma=$ARGV[2];
$n_box=$ARGV[3];
$input_file=$ARGV[4];
$out_file=$ARGV[5];


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
	$model[$i]=$data[0];
	for ($j=0;$j<9;$j++) {
	    $a[$i][$j][$n]=$data[1+2*$j];
	    $ia[$i][$j][$n]=$data[1+2*$j+1];
	    if (($j==$n_param)&&($i==$n_model)) {
		$a_fix[$n]=$a[$i][$j][$n];
		$s_y_out2[$n]=$a_fix[$n];
		if ($y_min>$a_fix[$n]) {
		    $y_min=$a_fix[$n];
		}
		if ($y_max<$a_fix[$n]) {
		    $y_max=$a_fix[$n];
		}
	    }
	}
    }
    $x[$n]=$n;
    $n++;
};
close(INPUT);


@s_y_out2=median_filter($n_box,\@a_fix);






pgbeg(0,'?',1,1);
pgscf(2.0);
pgenv($x[0],$x[$n-1],$y_min,$y_max,0,0);
pgsci(7);
pgline($n,\@x,\@a_fix);
pgsci(1);
pgpoint($n,\@x,\@a_fix,2);
#pgsci(15);
#pgline($n,\@x,\@s_y_out);
pgsci(2);
pgline($n,\@x,\@s_y_out2);
pgsci(1);
pgclos();
pgend();


open(OUTPUT,">$out_file");
for ($k=0;$k<$n;$k++) {
    print OUTPUT "$n_mod\n";
    for ($i=0;$i<$n_mod;$i++) {
	print OUTPUT "$model[$i] ";
	for ($j=0;$j<9;$j++) {
	    if (($j==$n_param)&&($i==$n_model)) {
		$a[$i][$j][$k]=$s_y_out2[$k];
	    }
	    print OUTPUT "$a[$i][$j][$k] $ia[$i][$j][$k] ";
	}
	print OUTPUT "\n";
    }
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



