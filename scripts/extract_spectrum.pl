#!/usr/bin/perl
#
#
#
use PGPLOT;


# This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s the and output file
# like the one created by "fit_3D_fg.pl"
# It extract a certain spectrum.
#
if ($#ARGV<6) {
    print "USE: extract_spec.pl MODEL N_MODEL START_WAVE STEP_WAVE INPUT_FILE SPEC_FILE N_SIGMA\n";        
    exit;
}

$model=$ARGV[0];
$n_model=$ARGV[1];
$start_w=$ARGV[2];
$step_w=$ARGV[3];
$input_file=$ARGV[4];
$spec_file=$ARGV[5];
$n_sigma=$ARGV[6];

#print "$n_sigma $#ARGV\n";
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
	if ($i==$n_model) {
	    if ($model eq "gauss") {
		$flux[$n]=2*3.1416*$a[$i][2][$n]*($a[$i][3][$n]**2);
#		$flux[$n]=2*$a[$i][2][$n]*$a[$i][3][$n];
#		$flux[$n]=$a[$i][2][$n];#*$a[$i][3][$n];
	    }
	    if ($model eq "disk") {
		$flux[$n]=2*3.1416*$a[$i][2][$n]*($a[$i][3][$n]**2);
#		$flux[$n]=2*$a[$i][2][$n]*$a[$i][3][$n];
#		$flux[$n]=$a[$i][2][$n];#*$a[$i][3][$n];
	    }
	    if ($model eq "vauc") {
		$flux[$n]=161280*2*3.1416*$a[$i][2][$n]*($a[$i][3][$n]**2);
#		$flux[$n]=2*$a[$i][2][$n]*$a[$i][3][$n];
#		$flux[$n]=$a[$i][2][$n];#*$a[$i][3][$n];
	    }
	}
    }
    if ($y_min>$flux[$n]) {
	$y_min=$flux[$n];
    }
    if ($y_max<$flux[$n]) {
	$y_max=$flux[$n];
    }
#    $w[$n]=$start_w+$step_w*($n+1);
    $w[$n]=$start_w+$step_w*($n);
    $n++;
};
close(INPUT);

$dev=$ENV{"PGPLOT_DEV"};

pgbeg(0,$dev,1,1);
pgscf(2.0);
pgenv($w[0],$w[$n-1],$y_min-0.05*($y_max-$y_min),$y_max+0.05*($y_max-$y_min),0,0);
pgline($n,\@w,\@flux);



for ($k=2;$k<$n-3;$k++) {
#    my @box;
    $j=0;
    for ($i=$k-2;$i<=$k+2;$i++) {
	if ($i!=$k) {
	    $box[$j]=$flux[$i];
#	    print "$j $box[$j] $i $flux[$i]\n";
	    $j++;
	}	
    }
#    print "@box = ";
    ($mean,$sigma)=stats(@box);
#    print "$mean+-$sigma\n";
    if (abs($flux[$k]-$mean)>($n_sigma*$sigma)) {
	$flux[$k]=$mean;
    }
};
pgsci(2);
pgline($n,\@w,\@flux);
pgsci(1);
pgclos();
pgend();


open(OUTPUT,">$spec_file");
#print OUTPUT "# $n 1 $start_w $w[$n-1] $step_w\n";
for ($k=0;$k<$n;$k++) {
    print OUTPUT "$k $w[$k] $flux[$k]\n";
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



