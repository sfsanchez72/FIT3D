#!/usr/bin/perl
#
#
#
use PGPLOT;


use PDL;

use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;

# This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s the and output file
# like the one created by "fit_3D_fg.pl"
# and fix one parameter to do
# a force fitting latter
if ($#ARGV<5) {
    print "USE: compare.pl PARAM N_MOD1 N_MOD2 N_MIN N_MAX INPUT_FILE\n";        
    exit;
}

$n_param=$ARGV[0];
$n_mod1=$ARGV[1];
$n_mod2=$ARGV[2];
$n_min=$ARGV[3];
$n_max=$ARGV[4];
$input_file=$ARGV[5];



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
	    if (($j==$n_param)&&($i==$n_mod1)) {
		$a_fix1[$n]=$a[$i][$j][$n];
	    }
	    if (($j==$n_param)&&($i==$n_mod2)) {
		$delta[$n]=$a[$i][$j][$n]-$a_fix1[$n];
		if ($y_min>$delta[$n]) {
		    $y_min=$delta[$n];
		}
		if ($y_max<$delta[$n]) {
		    $y_max=$delta[$n];
		}
	    }
	}
    }
    $x[$n]=$n;
    $n++;
};
close(INPUT);



# We fit to a polynomical function
# of the order npoly (using PDL);

($mean,$sigma)=stats(@delta);


pgbeg(0,'/XS',1,1);
pgscf(2.0);
pgenv($x[0],$x[$n-1],$y_min,$y_max,0,0);
pgpoint($n,\@x,\@delta,2);
pgclos();
pgend();

print "MEAN=$mean+-$sigma\n";

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


