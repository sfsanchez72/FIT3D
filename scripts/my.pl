#########################################
# My Subroutines
# S.F.Sanchez 2004
##########################################
$PDL::BIGPDL = 1;


#$naxes=read_naxes($fitsfile);
#@in_data=read_img($fitsfile);
#$npix=@$naxes[0];
#$nb_spec=@$naxes[1];
#write_img($outfile,$npix,$nb_spec,\@out_data); 

use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Primitive;
use lib qw(/z/ata/ssa/perl-5.8.0/lib /home/ssa/perl-5.8.0/lib );
use Socket;
#use PDL::GSLSF::LEGENDRE;

sub LINFIT1D {
   my $opthash = ref($_[-1]) eq "HASH" ? pop(@_) : {} ; 
   my %opt = parse( { Weights=>ones(1) }, $opthash ) ;
   barf "Usage: LINFIT1D incorrect args\n" if $#_<1 or $#_ > 3;
   my ($y, $fitfuncs, $wt) = @_;
   if ($#_ == 1) {
      ($y, $fitfuncs) = @_;
      $x = $y->xvals;
   }
   
#   my $wt = $opt{Weights};
   
   # Internally normalise data
   
   my $ymean = (abs($y)->sum)/($y->nelem);
   $ymean = 1 if $ymean == 0;
   my $y2 = $y / $ymean;
   
   # Do the fit
      
   my $M = $fitfuncs->xchg(0,1);
   my $C = $M->xchg(0,1) x ($M * $wt->dummy(0)) ;
   my $Y = $M->xchg(0,1) x ($y2->dummy(0) * $wt->dummy(0));

   # Fitted coefficients vector

   my   $a = matinv($C) x $Y;
   
   # Fitted data

   my $yfit = ($M x $a)->clump(2); # Remove first dim=1
   
   $yfit *= $ymean; # Un-normalise
   if (wantarray) {
      my $coeff = $a->clump(2);
      $coeff *= $ymean; # Un-normalise
      return ($yfit, $coeff);
   }
   else{
      return $yfit;
   }  
   
}
*LINFIT1D = \&LINFIT1D;

sub fit_legendre {
   my ($x, $y, $order) = @_;
   
   # means for each 1D data set
#   my $xmean = (abs($x)->average)->dummy(0);  # dummy for correct threading
#   my $ymean = (abs($y)->average)->dummy(0);
#   (my $tmp = $ymean->where($ymean == 0)) .= 1 if any $ymean == 0;
#   ($tmp = $xmean->where($xmean == 0)) .= 1 if any $xmean == 0;
   $xmean=1;
   $ymean=1;
   my $y2 = $y / $ymean;
   my $x2 = $x / $xmean;
   
   # Do the fit
      
   my $pow = sequence($order);
#   print $pow;
#   my $M = $x2->dummy(0) ** $pow;
   my @NX = $x2->dims;
   my $M = zeroes($order+1,$NX[0]);
   my $i,$j;
   for ($i=0;$i<$NX[0];$i++) {
       for ($j=0;$j<$order+1;$j++) {
	   my $valx=$x2->at($i);
	   my $val=legendre($j,$x2->at($i));
	   set($M,$j,$i,$val);
       }
   }


   my $C = $M->xchg(0,1) x ($M ) ;
   my $Y = $M->xchg(0,1) x ($y2->dummy(0) );

   # Fitted coefficients vector

   $a = matinv($C) x $Y;
   
   # Fitted data

   $yfit = ($M x $a)->clump(2); # Remove first dim=1
   
#   $yfit *= $ymean; # Un-normalise

   my $coeff = $a->clump(2);
#   for ($j=0;$j<$order+1;$j++) {
#       my $val = $coeff->at($j);
#       $YMEAN=$ymean->at(0);
#       $XMEAN=$xmean->at(0);
#       my $out = $YMEAN/legendre($j,$XMEAN);
#       $val=$val*$out;
#       set($coeff,$j,$val);
#   }
#   print "$yfit\n";

   return ($yfit, $coeff);
   
}




sub legendre {
    my $n=$_[0];
    my $x=$_[1];
    my $Pnx;
    my $Po;
    $Po=1;#legendre_norm($n,1000);
    $Pnx=legendre_norm($n,$x);
#    $Pnx=gsl_sf_legendre_Pl($x,$n);
    $Pnx=$Pnx/$Po;
    return $Pnx;
}

sub legendre_norm {
    my $n=$_[0];
    my $x=$_[1];
    my $Pnx;
    if ($n==0) {
	$Pnx=1;
    } else {
	if ($n==1) {
	    $Pnx=$x;
	} else {
	    $Pnx=((2*($n-1)+1)*$x*legendre_norm($n-1,$x)-($n-1)*legendre_norm($n-2,$x))/$n;
	}
	
    }
#    $Pnx=apr_n($Pnx,6);
    return $Pnx;
}


sub legendre_OLD {
    my $n=$_[0];
    my $x=$_[1];
    my $m;
    if ($n==(2*int($n/2))) {
	$m=$n/2;
    } else {
	$m=($n-1)/2;
    }
    my $Pnx=0;
    my $i;
    for ($i=0;$i<$m;$i++) {
	$Pnx=$Pnx+(((-1)**$i)*(fact(2*$n-2*$i))/((2**$n)*(fact($n-$i))*(fact($n-2*$i))))*($x**($n-2*$i));
    }
    return $Pnx;
}

sub fact {
    my $x=$_[0];
    my $f=1;
    for ($i=0;$i<$x;$i++) {
	$f=$f*($x-$i);
    }
    return $f;
}


sub gammaq {
    my $a=$_[0];
    my $x=$_[1];
    my $gamser,$gammcf,$gln;
    my $out;
   # print "$a $x\n";
    if (($x<0)||($a<=0)) {
	$out=0;
    } else {
	if ($x<($a+1)) {
	    $gamser=gser($a,$x);
	    $out=1-$gamser;
#	    print "GAMSER=$gamser OUT=$out\n";

	} else {
	    $gammcf=gcf($a,$x);
	 #   print "GAMMCF=$gammcf\n";
	    $out=$gammcf;
	}


    }
    return $out;
}

sub gser {
    my $a=$_[0];
    my $x=$_[1];
    my $ITMAX=100;
    my $EPS=3e-7;
    my $gln;
    my $n;
    my $sum,$del,$ap;
    $gln=gammaln($a);
    my $out;
#    print "X=$x,A=$a, gln=$gln\n";
    if ($x<=0) {
	$out=0;
    } else {
	$ap=$a;
	$sum=1.0/$a;
	$del=$sum;
	for ($n=1;$n<=$ITMAX;$n++) {
	    $ap=$ap+1;
	    $del=$del*$x/$ap;
	    $sum=$sum+$del;
	    my $AA=abs($del);
	    my $BB=$EPS*abs($sum);
#	    print "$n $AA<$BB\n";
	    if ($AA<$BB) {
		$out=$sum*exp(-$x+$a*log($x)-$gln);
#	    } else {
#		$n=$ITMAX+1;
	    }
#	    print "$out\n";
	}
    }
    return $out;
}


sub gcf {
    my $a=$_[0];
    my $x=$_[1];
    my $ITMAX=100;
    my $EPS=3e-7;
    my $FPMIN=1e-30;
    my $i,$an,$b,$c,$d,$del,$h;
    my $gln;
    $gln=gammaln($a);
    my $out;
    $b=$x+1-$a;
    $c=1/$FPMIN;
    $d=1/$b;
    $h=$d;
    for ($i=1;$i<=$ITMAX;$i++) {
	$an=-$i*($i-$a);
	$b=$b+2.0;
	$d=$an*$d+$b;
	if (abs($d)<$FPMIN) {
	    $d=$FPMIN;
	}
	$c=$b+$an/$c;
	if (abs($c)<$FPMIN) {
	    $c=$FPMIN;
	}
	$d=1/$d;
	$del=$d*$c;
	$h=$h*$del;
	if (abs($del-1.0)<$EPS) {
	    $i=$ITMAX+1;
	}
    }
    $out=exp(-$x+$a*log($x)-$gln)*$h;
    #print "OUT = $out $x $a $gln $h\n";
    return $out;
}




sub gamma {
    my $xx=$_[0];
    my $out=exp(gammaln($xx));
    return $out;
}

sub gammaln {
    my $xx=$_[0];
    my $x,$y,$tmp,$ser;
    my $out;
    my @cof=(76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5);
    my $j;
    if ($xx>0) {
	$x=$xx;
	$y=$x;
	$tmp=$x+5.5;
	$tmp=$tmp-($x+0.5)*log($tmp);
	$ser=1.000000000190015;
	for ($j=0;$j<6;$j++) {
	    $y=$y+1;
	    $ser=$ser+$cof[$j]/($y);
#	    print "SER=$ser $cof[$j] $y $j\n";
	}
	$out=-$tmp+log(2.5066282746310005*$ser/$x);
	$C=(2.5066282746310005*$ser/$x);
	$A=log(2.5066282746310005*$ser/$x);
#	print "OUT=$out,$tmp,$ser,$x,$A,$C\n";
    } else {
	$out=0;
    }
    return $out;

}


sub pmas_rss2maps {
    my $infile=$_[0];
    my $pre=$_[1];
    my $na=read_naxes($infile);    
    ($npix,$nb_spec)=@$na;
    my $n=sqrt($nb_spec);
    my $i,$j,$k;
    my @input_a=read_img($infile);
    for ($i=0;$i<$npix;$i++) {
	my @out_a;
	if ($i<10000) {
	    $text_i=$i;
	}
	if ($i<1000) {
	    $text_i="0".$i;
	}
	if ($i<100) {
	    $text_i="00".$i;
	}
	if ($i<10) {
	    $text_i="000".$i;
	}
	my $name=$pre."_".$text_i.".fits";
	for ($j=0;$j<$n;$j++) {
	    for ($k=0;$k<$n;$k++) {
		$out_a[$k][$j]=$input_a[$k+16*$j][$i];
	    }
	}
	system("rm $name");
	write_img($name,$n,$n,\@out_a);
    }
}


sub imcombine {
    my $w=$_[0];
    my $outfile=$_[1];
    my @images=@$w;
    my $j,$i,$n;
    my @array;
    my $crval1=read_img_header($images[0],"CRVAL1");
    my $cdelt1=read_img_header($images[0],"CDELT1");
#    $cd1_1=$read_img_header($images[0],"CD1_1");
#    $cd1_2=$read_img_header($images[0],"CD1_2");
    my $na=read_naxes($images[0]);    
    my $npix,$nb_spec;
    ($npix,$nb_spec)=@$na;
    my @a_tmp;
    for ($i=0;$i<($#images+1);$i++) {
	my @aa_tmp=read_img($images[$i]);
	for ($j=0;$j<$nb_spec;$j++) {
	    for ($k=0;$k<$npix;$k++) {
		for ($i=0;$i<($#images+1);$i++) {
		    $a_tmp[$j][$k][$i]=$aa_tmp[$j][$k];
		}
	    }
	}
    }

    for ($j=0;$j<$nb_spec;$j++) {
	for ($k=0;$k<$npix;$k++) {
	    my @f_tmp;
	    for ($i=0;$i<($#images+1);$i++) {
		$f_tmp[$i]=$a_tmp[$j][$k][$i];
#		print "$i $a_tmp[$j][$k][$i]";
	    }
	    $array[$j][$k]=median(@f_tmp);
#	    print "$j $k $array[$j][$k]\n";
	}
    }


#    for ($j=0;$j<$npix;$j++) {
#	for ($k=0;$k<$npix;$k++) {
#	    $array[$j][$k]/=$#images;
#	}
#    }
    system("rm $outfile");
    write_rss($outfile,$npix,$nb_spec,$crval1,$cdelt1,\@array);

}


sub extinction {
    my $w=$_[0];
    my @wa=@$w;
    my $n=$#wa+1;
    my $E_BV=$_[1];
#    my @E_L=(1000,1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000);
#    my @E_LV=(10.0,9.1,8.3,7.7,7.1,6.7,6.2,5.6,5.0,4.5,4.2,3.8,3.6,3.3,2.9,2.5);
    my $x=pdl(40000,26500,12200,6000,5470,4670,4110,2700,2600);
    my $y=pdl(0.000,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591);
    my $yi = interpol(pdl(@wa), $x, $y);
    my $yi = $yi;
    my $j,@out;
    for ($j=0;$j<$n;$j++) {
	$out[$j] = $E_BV*($yi->slice($j)->sclr);
	print "$j $wa[$j] $out[$j]\n";
    }
    return @out;
}

sub median_smooth {
  my $dat=$_[0];
  my $nx=$_[1];
  my $ny=$_[2];
  my $bx=$_[3];
  my $by=$_[4];
  my $i,$j,$ii,$jj;
  my @data=@$dat;
  my @out_data;
#  print "NN=$nx,$ny $bx,$by\n";  
#  for ($i=$bx;$i<$nx-$bx;$i++) {
#      for ($j=$by;$j<$ny-$by;$j++) {
  for ($i=0;$i<$nx;$i++) {
      for ($j=0;$j<$ny;$j++) {
	  $out_data[$j][$i]=$data[$j][$i];
	  if (($i>$bx)&&($i<$nx-$bx)&&($j>$by)&&($j<$ny-$by)) {
	      my @my_tmp;
	      my $nn=0;
	      for ($ii=-$bx;$ii<$bx;$ii++) {
		  for ($jj=-$by;$jj<$by;$jj++) {
		      $my_tmp[$nn]=$data[$j+$jj][$i+$ii];
		      $nn++;
		  }
	      }
	      $out_data[$j][$i]=median(@my_tmp);
#	      print "$j $i $data[$j][$i] $out_data[$j][$i]\n";
	  }
      }
# sort numerically ascending
  }
  return @out_data;
}


sub mean_smooth {
  my $dat=$_[0];
  my $nx=$_[1];
  my $ny=$_[2];
  my $bx=$_[3];
  my $by=$_[4];
  my $i,$j,$ii,$jj;
  my @data=@$dat;
  my @out_data;
#  print "NN=$nx,$ny $bx,$by\n";  
#  for ($i=$bx;$i<$nx-$bx;$i++) {
#      for ($j=$by;$j<$ny-$by;$j++) {
  for ($i=0;$i<$nx;$i++) {
      for ($j=0;$j<$ny;$j++) {
	  $out_data[$j][$i]=$data[$j][$i];
	  if (($i>$bx)&&($i<$nx-$bx)&&($j>$by)&&($j<$ny-$by)) {
	      my @my_tmp;
	      my $nn=0;
	      for ($ii=-$bx;$ii<$bx;$ii++) {
		  for ($jj=-$by;$jj<$by;$jj++) {
		      $my_tmp[$nn]=$data[$j+$jj][$i+$ii];
		      $nn++;
		  }
	      }
	      $out_data[$j][$i]=mean(@my_tmp);
#	      print "$j $i $data[$j][$i] $out_data[$j][$i]\n";
	  }
      }
# sort numerically ascending
  }
  return @out_data;
}


sub median_clip {
  my $dat=$_[0];
  my $nx=$_[1];
  my $ny=$_[2];
  my $bx=$_[3];
  my $by=$_[4];
  my $nsig=$_[5];
  my $i,$j,$ii,$jj;
  my @data=@$dat;
  my @out_data;
#  for ($i=$bx;$i<$nx-$bx;$i++) {
#      for ($j=$by;$j<$ny-$by;$j++) {
  for ($i=0;$i<$nx;$i++) {
      for ($j=0;$j<$ny;$j++) {
	  $out_data[$j][$i]=$data[$j][$i];
	  if (($i>$bx)&&($i<$n1-$bx)&&($j>$by)&&($j<$n2-$by)) {
	      my @my_tmp;
	      my $nn=0;
	      for ($ii=-$bx;$ii<$bx;$ii++) {
		  for ($jj=-$by;$jj<$by;$jj++) {
		      $my_tmp[$nn]=$data[$j+$jj][$i+$ii];
		      $nn++;
		  }
	      }
#	      $out_data[$j][$i]=median(@my_tmp);
	      my $tmp_med=median(@my_tmp);
	      my $tmp_sig=sigma(@my_tmp);
	      if (abs($data[$j][$i]-$tmp_med)>($nsig*$tmp_sig)) {
		  $out_data[$j][$i]=$tmp_med;
	      }
#	      print "$j $i $data[$j][$i] $out_data[$j][$i]\n";
	  }
      }
# sort numerically ascending
  }
  return @out_data;
}

sub median_clip1d {
  my $dat=$_[0];
  my $nx=$_[1];
  my $bx=$_[2];
  my $nsig=$_[3];
  my $i,$j,$ii,$jj;
  my @data=@$dat;
  my @out_data;
  for ($i=0;$i<$nx;$i++) {
      $out_data[$i]=$data[$i];
      if (($i>$bx)&&($i<$nx-$bx)) {
	  my @my_tmp;
	  my $nn=0;
	  for ($ii=-$bx;$ii<$bx;$ii++) {
	      $my_tmp[$nn]=$data[$i+$ii];
	      $nn++;
	  }
	  my $tmp_med=median(@my_tmp);
	  my $tmp_sig=sigma(@my_tmp);
#	  print "MY=$i $data[$i] $out_data[$i] $tmp_med $tmp_sig $nsig\n";

	  if (abs($data[$i]-$tmp_med)>($nsig*$tmp_sig)) {
	      $out_data[$i]=$tmp_med;
	  }
      }

  }
  return @out_data;
}




sub median {
  my @data=@_;
  my @out_data = sort {$a <=> $b} @data;
  my $N=$#out_data+1;
  my $median_sub;
  if ($N==2*(int($N/2))) {
      $median_sub=0.5*($out_data[$N/2]+$out_data[$N/2+1]);
  } else {
      $median_sub=$out_data[($N+1)/2];
  }
#  print "N=$N $median\n";
  return $median_sub;
}

sub median_cut {
  my $in_data=$_[0];
  my $min=$_[1];
  my $max=$_[2];
  my @in_array=@$in_data;
  my $k=0;
  for ($j=0;$j<($#in_array+1);$j++) {
      if (($in_array[$j]>$min)&&($in_array[$j]<$max)) {
	  $data[$k]=$in_array[$j];
	  $k++;
	  }
  }
# sort numerically ascending
  my @out_data = sort {$a <=> $b} @data;
  my $n_median=int(($#out_data+1)*0.5);
  my $median = ($out_data[$n_median-1]+$out_data[$n_median]+$out_data[$n_median+1])/3;
  return $median;
}

sub median_filter {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    my @val=@$array;
    my $i,$jk;
    for ($i=$box;$i<($#val+1-$box);$i++) {
	my @tmp;
	for ($jk=0;$jk<2*$box;$jk++) {
	    $tmp[$jk]=@$array[$i-$box+$jk];
	}
	$val[$i]=median(@tmp);
    }

    for ($i=1;$i<$box;$i++) {
	my @tmp;
	$effec_box=$i;
	for ($jk=0;$jk<2*$effec_box;$jk++) {
	    $tmp[$jk]=@$array[$i-$effec_box+$jk];
	}
	$val[$i]=median(@tmp);
    }

    for ($i=($#val+1-$box);$i<$#val;$i++) {
	my @tmp;
	$effec_box=$#val-$i+1;
	for ($jk=0;$jk<2*$effec_box;$jk++) {
	    $tmp[$jk]=@$array[$i-$effec_box+$jk];
	}
	$val[$i]=median(@tmp);
#	$val[$i]=$val[$#val-$box-1];
    }
    return @val;
}


sub smotth_filter {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    my @val=@$array;
    my $i,$jk;
    for ($i=$box;$i<($#val+1-$box);$i++) {
	my @tmp;
	for ($jk=0;$jk<2*$box;$jk++) {
	    $tmp[$jk]=@$array[$i-$box+$jk];
	}
	$val[$i]=median(@tmp);
    }

    for ($i=1;$i<$box;$i++) {
	my @tmp;
	$effec_box=$i;
	for ($jk=0;$jk<2*$effec_box;$jk++) {
	    $tmp[$jk]=@$array[$i-$effec_box+$jk];
	}
	$val[$i]=median(@tmp);
    }

    for ($i=($#val+1-$box);$i<$#val;$i++) {
	my @tmp;
	$effec_box=$#val-$i+1;
	for ($jk=0;$jk<2*$effec_box;$jk++) {
	    $tmp[$jk]=@$array[$i-$effec_box+$jk];
	}
	$val[$i]=median(@tmp);
#	$val[$i]=$val[$#val-$box-1];
    }
    return @val;
}


sub median_filterB {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    my @val=@$array;
    my $i,$jk;
    for ($i=$box;$i<($#val+1-$box);$i++) {
	my @tmp;
	for ($jk=0;$jk<2*$box;$jk++) {
	    $tmp[$jk]=@$array[$i-$box+$jk];
	}
	$val[$i]=median(@tmp);
    }

    return @val;
}


sub median_box {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    my @in_val=@$array;
    my @out_val;
    my $i,$j;
    my $k=0;
    for ($i=$box;$i<($#in_val+1-$box);$i=$i+2*$box) {
	my @tmp;
	for ($j=0;$j<2*$box;$j++) {
	    $tmp[$j]=@$array[$i-$box+$j];
	}
	$out_val[$k]=median(@tmp);
	$k++;
    }
    return @out_val;
}

sub mean_box {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    @in_val=@$array;
    my $i,$j,$k;
    $k=0;
    for ($i=$box;$i<($#in_val+1-$box);$i=$i+2*$box) {
	my @tmp;
	for ($j=0;$j<2*$box;$j++) {
	    $tmp[$j]=@$array[$i-$box+$j];
	}
	$out_val[$k]=mean(@tmp);
#	$out_sval[$k]=sigma(@tmp);
	$k++;
    }
    return @out_val;
}


sub minmax {
  my @data=@_;
  my $sum=0; 
  my $j;
  my $min=1e20;
  my $max=-1e20;
  for ($j=0;$j<($#data+1);$j++) {
      if ($min>$data[$j]) {
	  $min=$data[$j];	 
      }
      if ($max<$data[$j]) {
	  $max=$data[$j];	  
      }
  }
  return ($min,$max);
}

sub mean {
  my @data=@_;
  my $sum=0; 
  my $j;
  my $mean;
  for ($j=0;$j<($#data+1);$j++) {
      $sum=$sum+$data[$j];
#      print "MEAN = $j $sum $data[$j]\n";
  }
  if ($#data>-1) {
      $mean = $sum/($#data+1);
  } else {
      $mean = 0;
  }
  return $mean;
}

sub sigma {
  my @data=@_;
  my $mean = mean(@data);
  my $stddev = 0;
  my $j;
  my $sum=0; 
  for ($j=0;$j<($#data+1);$j++) {
#      print "$j $sum $data[$j] $mean\n";
      $sum=$sum+($data[$j]-$mean)**2;
  }
  if ($#data>0) {
      $sum=$sum/($#data);
      $stddev=sqrt($sum);
  } else {
      $stddev=0;
  }
#  print "SIGMA = $stddev $mean $#data\n";

  return $stddev;
}


sub sigma_abs {
  my @data=@_;
  my $mean = mean(@data);
  my $stddev = 0;
  my $j;
  my $sum=0;
  my @adata;
  for ($j=0;$j<($#data+1);$j++) {
      $adata[$j]=abs($data[$j]-$mean);      
  }
  $stddev=median(@adata);
  return $stddev;
}


sub sigma_m {
  my @data=@_;
  my $mean = median(@data);
  my $stddev = 0;
  my $j;
  my $sum=0; 
  for ($j=0;$j<($#data+1);$j++) {
      $sum=$sum+($data[$j]-$mean)**2;
  }
  if ($#data>0) {
      $sum=$sum/($#data);
      $stddev=sqrt($sum);
  } else {
      $stddev=0;
  }
  return $stddev;
}

sub sigma_m_pos {
  my @data=@_;

  my @new_data;
  my $i;
  my $k=0;
  for ($i=0;$i<($#data+1);$i++) {
      if ($data[$i]>0) {
	  $new_data[$k]=$data[$i];
	  $k++;
      }
  }


  my $mean = median(@new_data);
  my $stddev = 0;
  my $j;
  my $sum=0; 
  for ($j=0;$j<($#new_data+1);$j++) {
      $sum=$sum+($new_data[$j]-$mean)**2;
  }
  if ($#new_data>0) {
      $sum=$sum/($#new_data);
      $stddev=sqrt($sum);
  } else {
      $stddev=0;
  }
  return $stddev;
}


sub read_naxes() {
    my $file=@_[0];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    $fptr->close_file($status);
    return $naxes;
}

sub read_img_header() {
    my $file=$_[0];
    my $header_id=$_[1];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    my $value,$comment;
    $fptr->read_keyword($header_id,$value,$comment,$status);
    $fptr->close_file($status);
    return $value;
}

sub write_img_header() {
    my $file=$_[0];
    my $header_id=$_[1];
    my $value=$_[2];
    my $type=$_[3];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,READWRITE,$status);
    if ($type eq "TBYTE") {
	$fptr->write_key_byt($header_id,$value,$header_id,$status);
    }
    if ($type eq "TSHORT") {
	$fptr->write_key_sht($header_id,$value,$header_id,$status);
    }
    if ($type eq "TLONG") {
	$fptr->write_key_lng($header_id,$value,$header_id,$status);
    }
    if ($type eq "TFLOAT") {
	$fptr->write_key_flt($header_id,$value,$header_id,$status);
    }
    if ($type eq "TDOUBLE") {
	$fptr->write_key_dbl($header_id,$value,$header_id,$status);
    }
    if ($type eq "TSTRING") {
	$fptr->write_key_str($header_id,$value,$header_id,$status);
    }
    $fptr->close_file($status);
    return $value;
}

sub read_img_headers() {
    my $file=$_[0];
    my $tmp=$_[1];
    my @headers_id=@$tmp;
#    print "$file @headers_id\n";
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    my $value,$comment;
    my @values;
    my $i;
    for ($i=0;$i<$#headers_id+1;$i++) {
	my $header_id=$headers_id[$i];
#	print "'$header_id'\n";
	$fptr->read_keyword($header_id,$value,$comment,$status);
	$values[$i]=$value;
#	print "'$header_id' $values[$i] $value\n";
    }
    $fptr->close_file($status);
    return @values;
}



sub read_img() {
    my $file=@_[0];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    my $naxis1;
    my $naxis2;
    ($naxis1,$naxis2) = @$naxes;
    my (@array, $nullarray, $anynull);
    $fptr->read_pixnull(Astro::FITS::CFITSIO::TFLOAT(), [1,1], $naxis1*$naxis2, \@array, $nullarray, $anynull ,$status);
    $fptr->close_file($status);
    return @array;
}

sub read_img_sec() {
    my $file=@_[0];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    my $naxis1;
    my $naxis2;
    ($naxis1,$naxis2) = @$naxes;
    my (@array, $nullarray, $anynull);
    $fptr->read_pixnull(Astro::FITS::CFITSIO::TFLOAT(), [1,1], $naxis1*$naxis2, \@array, $nullarray, $anynull ,$status);
    $fptr->close_file($status);
    return @array;
}


sub read_img_int() {
    my $file=@_[0];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    my $naxis1;
    my $naxis2;
    ($naxis1,$naxis2) = @$naxes;
    my (@array, $nullarray, $anynull);
    $fptr->read_pixnull(Astro::FITS::CFITSIO::TINT(), [1,1], $naxis1*$naxis2, \@array, $nullarray, $anynull ,$status);
    $fptr->close_file($status);
    return @array;
}


sub read_line() {
    my $file=@_[0];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    my $naxis1;
    ($naxis1) = @$naxes;
    my (@array, $nullarray, $anynull);
    $fptr->read_img(TFLOAT,1,$naxis1,$nullarray,\@array,$anynull,$status);
#    $fptr->read_pixnull(Astro::FITS::CFITSIO::TFLOAT(), [1], $naxis1, \@array, $nullarray, $anynull ,$status);
    $fptr->close_file($status);
    return @array;
}



sub read_cube() {
    my $file=@_[0];
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    $fptr->get_img_parm(undef,undef,$naxes,$status);
    my $naxis1;
    my $naxis2;
    my $naxis3;
    ($naxis1,$naxis2,$naxis3) = @$naxes;
    my (@array, $nullarray, $anynull);
    $fptr->read_pixnull(Astro::FITS::CFITSIO::TFLOAT(), [1,1,1], $naxis1*$naxis2*$naxis3, \@array, $nullarray, $anynull ,$status);
    $fptr->close_file($status);
    return @array;
}

sub write_img() {
     my $outfile=$_[0];
     my $npix=$_[1];
     my $nb_spec=$_[2];
     my $array=$_[3];
     my $simple=1;
     my $bitpix=-32;
     my $naxis=2;
     my $naxes2=[$npix,$nb_spec];
     my $status =0;
     @val=@$array;
     my $out=Astro::FITS::CFITSIO::create_file($outfile,$status);
     $out->write_imghdr($bitpix,$naxis,$naxes2,$status);
     $out->write_img(Astro::FITS::CFITSIO::TFLOAT(),1, $nb_spec*$npix+1, \@val, $status);
     $out->close_file($status);
}


sub update_fits() {
     my $outfile=$_[0];
     my $naxes2=$_[1];
     my $naxis=$_[2];
     my $array=$_[3];
     my $simple=1;
     my $bitpix=-32;
     my $status =0;
     my @val=@$array;
     @axes=@$naxes2;
     my $npix=1;
     my $i,$j,$k;
     for ($i=0;$i<$naxis;$i++) {
	 $npix=$npix*$axes[$i];
     }     
     my $out=Astro::FITS::CFITSIO::open_file($outfile,READWRITE,$status);
#     my $out=Astro::FITS::CFITSIO::create_file($outfile,$status);
#     $out->write_imghdr($bitpix,$naxis,$naxes2,$status);
     $out->write_img(Astro::FITS::CFITSIO::TFLOAT(),1, $npix+1, \@val, $status);     $out->close_file($status);
     $out->close_file($status);
}





#
# write_fits (including cubes)
# write_fits (out_file,[n_i],Num.axis,@array);
sub write_fits() {
     my $outfile=$_[0];
     my $naxes2=$_[1];
     my $naxis=$_[2];
     my $array=$_[3];
     my $simple=1;
     my $bitpix=-32;
     my $status =0;
     my @val=@$array;
     @axes=@$naxes2;
     my $npix=1;
     my $i,$j,$k;
     for ($i=0;$i<$naxis;$i++) {
	 $npix=$npix*$axes[$i];
     }     
     my $out=Astro::FITS::CFITSIO::create_file($outfile,$status);
     $out->write_imghdr($bitpix,$naxis,$naxes2,$status);
     $out->write_img(Astro::FITS::CFITSIO::TFLOAT(),1, $npix+1, \@val, $status);     $out->close_file($status);
     $out->close_file($status);
}




sub write_rss() {
     my $outfile=$_[0];
     my $npix=$_[1];
     my $nb_spec=$_[2];
     my $start=$_[3];
     my $delta=$_[4];
     my $array=$_[5];
     my $simple=1;
     my $bitpix=-32;
     my $naxis=2;
     my $naxes2=[$npix,$nb_spec];
     my $status =0;
     @val=@$array;
     my $out=Astro::FITS::CFITSIO::create_file($outfile,$status);
     $out->write_imghdr($bitpix,$naxis,$naxes2,$status);
     $out->write_key_lng('CRPIX1',1,'CRpix1',$status) and print "CRPIX1 status = $status\n";
     $out->write_key_lng('CRPIX2',1,'CRpix2',$status) and print "CRPIX1 status = $status\n";
     $out->write_key_flt('CRVAL1',$start,5,'CRVAL',$status) and print "CRVAL status = $status\n";
     $out->write_key_flt('CRVAL2',1,5,'CRVAL',$status) and print "CRVAL status = $status\n";
     $out->write_key_flt('CDELT1',$delta,5,'Cdelta1',$status) and print "CDELT1 status = $status\n";
     $out->write_key_flt('CDELT2',1,5,'Cdelta2',$status) and print "CDELT1 status = $status\n";
     $out->write_key_str('CTYPE1','WAVELENGTH','CTYPE1',$status) and print "CTYPE status = $status\n";
     $out->write_key_str('CTYPE2','SCAN','CTYPE2',$status) and print "CTYPE status = $status\n";
     $out->write_key_flt('CD1_1',$delta,5,'Cd1_1',$status) and print "CD1_1 status = $status\n";
     $out->write_key_flt('CD1_2',1,5,'Cd1_2',$status) and print "CD1_2 status = $status\n";
     $out->write_img(Astro::FITS::CFITSIO::TFLOAT(),1, $nb_spec*$npix+1, \@val, $status);
     $out->close_file($status);
}


sub write_wcs() {
     my $outfile=$_[0];
     my $crpix1=$_[1];
     my $crval1=$_[2];
     my $ctype1=$_[3];
     my $cd1_1=$_[4];
     my $cd2_1=$_[5];
     my $crpix2=$_[6];
     my $crval2=$_[7];
     my $ctype2=$_[8];
     my $cd1_2=$_[9];
     my $cd2_2=$_[10];
     my $simple=1;
     my $bitpix=-32;
     my $naxis=2;     
     my $status =0;     
     my $out=Astro::FITS::CFITSIO::open_file($outfile,READWRITE,$status);
     $out->write_key_lng('CRPIX1',$crpix1,'CRpix1',$status) and print "CRPIX1 status = $status\n";
     $out->write_key_flt('CRVAL1',$crval1,5,'CRVAL1',$status) and print "CRVAL status = $status\n";
     $out->write_key_str('CTYPE1',$ctype1,'CTYPE1',$status) and print "CTYPE status = $status\n";
     $out->write_key_flt('CD1_1',$cd1_1,5,'Cd1_1',$status) and print "CD1_1 status = $status\n";
     $out->write_key_flt('CD2_1',$cd2_1,5,'Cd2_1',$status) and print "CD1_2 status = $status\n";
     $out->write_key_lng('CRPIX2',$crpix2,'CRpix2',$status) and print "CRPIX1 status = $status\n";
     $out->write_key_flt('CRVAL2',$crval2,5,'CRVAL2',$status) and print "CRVAL status = $status\n";
     $out->write_key_str('CTYPE2',$ctype2,'CTYPE2',$status) and print "CTYPE status = $status\n";
     $out->write_key_flt('CD1_2',$cd1_2,5,'Cd1_2',$status) and print "CD1_1 status = $status\n";
     $out->write_key_flt('CD2_2',$cd2_2,5,'Cd2_2',$status) and print "CD1_2 status = $status\n";
     $out->close_file($status);
}

sub write_crval_cube() {
     my $outfile=$_[0];
     my $tmp=$_[1];
     my @data=@$tmp;
     my $simple=1;
     my $bitpix=-32;
     my $naxis=2;     
     my $status =0;     
     my $out=Astro::FITS::CFITSIO::open_file($outfile,READWRITE,$status);
     $out->write_key_lng('CRPIX1',$data[0],'CRpix1',$status) and print "CRPIX1 status = $status\n";
     $out->write_key_flt('CRVAL1',$data[1],5,'CRVAL1',$status) and print "CRVAL1 status = $status\n";
     $out->write_key_flt('CDELT1',$data[2],5,'CDELT1',$status) and print "CDELT1 status = $status\n";
     $out->write_key_lng('CRPIX2',$data[3],'CRpix2',$status) and print "CRPIX2 status = $status\n";
     $out->write_key_flt('CRVAL2',$data[4],5,'CRVAL2',$status) and print "CRVAL2 status = $status\n";
     $out->write_key_flt('CDELT2',$data[5],5,'CDELT2',$status) and print "CDELT2 status = $status\n";
     $out->write_key_lng('CRPIX3',$data[6],'CRpix3',$status) and print "CRPIX3 status = $status\n";
     $out->write_key_flt('CRVAL3',$data[7],5,'CRVAL3',$status) and print "CRVAL3 status = $status\n";
     $out->write_key_flt('CDELT3',$data[8],5,'CDELT3',$status) and print "CDELT3 status = $status\n";
     $out->close_file($status);
}





sub log10 {
  my $n = shift;
  if ($n==0) {
      $n=1e-12;
  }
  return log(abs($n))/log(10);
}

sub hist_y {
    local ($start,$end,$nbin,@data)=@_;
    my @out;
    my $i;
    my $j;
#    print "$n_s, $start, $end, $nbin\n";
    for ($i=0;$i<$nbin;$i++) {
	for ($j=0;$j<($#data+1);$j++) {
	    if ($j==0) {
		$out[$i]=0
		} else {
		    if (($data[$j]>=($start+(($end-$start)/$nbin)*$i))&&($data[$j]<($start+(($end-$start)/$nbin)*($i+1)))) {
			$out[$i]=$out[$i]+1;
		    }
		}
	}
#	$out[$i]=($out[$i]/$#data);
    }
    return @out;
}

sub hist_y_norm {
    local ($start,$end,$nbin,@data)=@_;
    my @out;
    my $i;
    my $j;
    my $count;
#    print "$n_s, $start, $end, $nbin\n";
    for ($i=0;$i<$nbin;$i++) {
	for ($j=0;$j<($#data+1);$j++) {
	    if ($j==0) {
		$out[$i]=0
		} else {
		    if (($data[$j]>=($start+(($end-$start)/$nbin)*$i))&&($data[$j]<($start+(($end-$start)/$nbin)*($i+1)))) {
			$out[$i]=$out[$i]+1;
			$count++;
		    }
		}
	}
#	$out[$i]=($out[$i]/$count);
    }
    for ($i=0;$i<$nbin;$i++) {
	if ($count>0) {
	    $out[$i]=($out[$i]/$count);
	} else {
	    $out[$i]=0;
	}
    }
    return @out;
}


sub hist_y_norm_s {
    local ($start,$end,$nbin,@data)=@_;
    my @out;
    my $i;
    my $j;
    my $count;
#    print "$n_s, $start, $end, $nbin\n";
    $out[0]=0;
    for ($i=1;$i<=$nbin;$i++) {
	for ($j=0;$j<($#data+1);$j++) {
	    if ($j==0) {
		$out[$i]=0
		} else {
		    if (($data[$j]>=($start+(($end-$start)/$nbin)*$i))&&($data[$j]<($start+(($end-$start)/$nbin)*($i+1)))) {
			$out[$i]=$out[$i]+1;
			$count++;
		    }
		}
	}
#	$out[$i]=($out[$i]/$count);
    }
    for ($i=1;$i<=$nbin;$i++) {
	if ($count>0) {
	    $out[$i]=($out[$i]/$count);
	} else {
	    $out[$i]=0;
	}
    }
    return @out;
}

sub hist_x_s {
    local ($start,$end,$nbin)=@_;
    my @out;
    my $i;
    my $j;
    for ($i=0;$i<=$nbin;$i++) {
	$out[$i]=$start+(($end-$start)/$nbin)*($i+0.5);			
    }
    return @out;
}

sub hist_x {
    local ($start,$end,$nbin)=@_;
    my @out;
    my $i;
    my $j;
    for ($i=0;$i<$nbin;$i++) {
	$out[$i]=$start+(($end-$start)/$nbin)*($i+0.5);			
    }
    return @out;
}



sub l_abs {
    my @data=@_;
    my $ho=$data[0];
    my $qo=$data[1];
    my $z=$data[2];
#    my $alf=$data[3];
    my $m=$data[3];
#    my $kcor=$data[4];
    my $dl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
#    my $m_gal=$m*(3*(10**10)*$dl)**2;
#    my $m_gal=$m*((1+$z)**(2))*4*3.14159*(3*(10**5)*$dl*3.09e24)**2;
#    my $m_gal=$m*((1+$z))*4*3.14159*(3*(10**5)*$dl*3.09e24)**2;
    my $m_gal=$m*4*3.14159*(3*(10**5)*$dl*3.09e24)**2;
#    my $m_gal=$m*4*3.14159*(3*(10**5)*$dl*3.09e24)**2;
    return $m_gal;
}

sub m_abs {
    my @data=@_;
    my $ho=$data[0];
    my $qo=$data[1];
    my $z=$data[2];
    my $m=$data[3];
    my $kcor=$data[4];
    my $dl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
    my $m_gal=$m+5-5*log10(3*(10**11)*$dl)-$kcor;
    return $m_gal;
}

sub m_abs_low {
    my @data=@_;
    my $ho=$data[0];
    my $z=$data[1];
    my $m=$data[2];
    my $dl=$z/$ho;
    my $m_gal=$m+5-5*log10(3*(10**9)*$dl);
    return $m_gal;
}

sub m_obs {
    my @data=@_;
    my $ho=$data[0];
    my $qo=$data[1];
    my $z=$data[2];
    my $m=$data[3];
    my $kcor=$data[4];
    my $dl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
    my $m_gal=$m-5+5*log10(3*(10**11)*$dl)+$kcor;
    return $m_gal;
}

sub re_abs {
    my @data=@_;
    my $ho=$data[0];
    my $qo=$data[1];
    my $z=$data[2];
    my $re=$data[3];
    my $kcor=$data[4];
    my $dl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
    my $da=$dl/((1+$z)**2);
    my $re_k=1454.44*$da*$re;
    return $re_k;
}

sub d_abs {
    my @data=@_;
    my $ho=$data[0];
    my $qo=$data[1];
    my $z=$data[2];
    my $re=$data[3];
    my $kcor=$data[4];
    my $dl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
   # my $da=$dl/((1+$z)**2);
    my $re_k=1454.44*$dl*$re;
    return $re_k;
}

sub apr {
    my $z=@_[0];
    my @s=split(/\./,$z);
    my $last=substr($s[1],0,-length($s[1])+2);
    my $zz=$s[0].".".$last;
    return $zz;
}

sub apr_n {
    my $z=$_[0];
    my $n=$_[1];
    my $zz;
    my @s=split(/\./,$z);
    my $last=substr($s[1],0,-length($s[1])+$n);
     if ($last) {
	$zz=$s[0].".".$last;
    } else {
	$zz=$s[0];
    }
    return $zz;
}

#
# Cosmology
#






sub invzdot {
    my $z = $_[0];
    my $om = $_[1];
    my $ol =$_[2];
    my $inv=(1+$z)*(1+$z)*sqrt(1+$om*$z+$ol*(-1+(1+$z)**(-2)));
    my $out=1/$inv;
    return $out;
}

sub tint {
    my $zmax=$_[0];
    my $om=$_[1];
    my $ol=$_[2];
    my $sum=0;
    my $zi=0;
    my $dz=0.001;
    my $inc;
    my $sum=0;
    while ($zi <= ($zmax-$dz)) {
	$inc = (invzdot($zi,$om,$ol)+invzdot($zi+$dz,$om,$ol))/2;
	$sum+=$dz*$inc;
	$zi+=$dz;
    }
    return $sum;
}

sub chi {
    my $z=$_[0];
    my $om=$_[1];
    my $ol=$_[2];
    my $kk;
    if (($om+$ol)==1) {
	$kk=1;
    } else  {
	$kk=sqrt(abs(1-$om-$ol));
    }
    my $ch=$kk/((1+$z)*sqrt($om*$z+1+$ol*((1+$z)**(-2)-1)));
    return $ch;
}

sub chiint {
    my $zmax=$_[0];
    my $om=$_[1];
    my $ol=$_[2];
    my $dz=0.001;
    my $sum=0;
    my $zi=0;
    while ($zi <= ($zmax-$dz)) {
	$inc = (chi($zi,$om,$ol) + chi($zi+$dz,$om,$ol))/2;
	$sum += $dz*$inc;
	$zi += $dz;
    }
    return $sum;
}

sub sinh2 {
    my $x=$_[0];
    my $y=(exp($x)-exp(-$x))/2;
    return $y;
}

$cosmo = {
     ludis => 0, # Luminosity Distance (Mpc)
     angsize => 0, # Angular Distance (Mpc)
     propmo => 0, # Proper Motion Distance (Mpc)
     lbtime => 0, # Look-back time (Gyr)
     dm => 0, # Distance Modulus
     scale => 0 # Angular Scale (Kpc/");
 };

sub cosm {
    my $z=$_[0];
    my $h=$_[1];
    my $om=$_[2];
    my $ol=$_[3];
    my $c=2.9979e5;
    my $th=978.0/$h;
    my $ch=chiint($z,$om,$ol);
    my $k,$kk,$sigma,$r0,$r;
    if (($om+$ol)==1) {
	$k=0; 
	$sigma=$ch;
	$r0=$c/$h;
    }
    if (($om+$ol)>1) {
	$k=1;
	$sigma=sin($ch);
	$kk=sqrt($om+$ol-1);
	$r0=$c/($h*$kk);
    }
    if (($om+$ol)<1) {
	$k=-1;
	$kk=sqrt(1-$om-$ol);
	$r0=$c/($h*$kk);
	$sigma=sinh2($ch);
    }

    $r=($r0/(1+$z))*$sigma;
    if ($r<1) {
	$r=1;
    }
    my $ludis=$r*((1+$z)**2);
    my $angsize=$r;
    my $propmo=$r*(1+$z);
    my $lbtime=$th*tint($z,$om,$ol);
    my $dm=5*0.4343*log($r*((1+$z)**2))+25;
    my $scale=$angsize/206.264806;

    $cosmo->{ludis}=$ludis;
    $cosmo->{angsize}=$angsize;
    $cosmo->{propmo}=$propmo;
    $cosmo->{lbtime}=$lbtime;
    $cosmo->{dm}=$dm;
    $cosmo->{scale}=$scale; # Angular Scale (Kpc/");
 
#    print "$scale\n";
    return ($ludis,$angsize,$propmo,$lbtime,$dm,$scale);
}


sub m_abs_lambda {
    my @data=@_;
    my $ho=$data[0];
    my $om=$data[1];
    my $ol=$data[2];
    my $z=$data[3];
    my $m=$data[4];
    my $kcor=$data[5];
    
    my @cosmo=cosm($z,$ho,$om,$ol);

#    my $qo=($om+$ol)/2;
#    my $dl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
#    my $odl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
    my $dl=$cosmo[0];
#    my $rat=$odl/$dl;
#    print "$odl|$dl|$rat\n";
   # print "$m+5-5*log10(3*(10**11)*$dl)-$kcor $qo\n";
#    my $m_gal=$m+5-5*log10(3*(10**11)*$dl)-$kcor;
    my $m_gal=$m-5*log10($dl*100000)-$kcor;
    return $m_gal;
}

sub m_obs_lambda {
    my @data=@_;
    my $ho=$data[0];
    my $om=$data[1];
    my $ol=$data[2];
    my $z=$data[3];
    my $m=$data[4];
    my $kcor=$data[5];
    
    my @cosmo=cosm($z,$ho,$om,$ol);

#    my $qo=($om+$ol)/2;
#    my $dl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
#    my $odl=(1/($ho*$qo**2))*(($z)*$qo+($qo-1)*((2*$qo*($z)+1)**0.5-1));
    my $dl=$cosmo[0];
#    my $rat=$odl/$dl;
#    print "$odl|$dl|$rat\n";
   # print "$m+5-5*log10(3*(10**11)*$dl)-$kcor $qo\n";
#    my $m_gal=$m+5-5*log10(3*(10**11)*$dl)-$kcor;
    my $m_gal=$m+5*log10($dl*100000)+$kcor;
    return $m_gal;
}


sub re_abs_lambda {
    my @data=@_;
    my $ho=$data[0];
    my $om=$data[1];
    my $ol=$data[2];
    my $z=$data[3];
    my $re=$data[4];
    my @cosmo=cosm($z,$ho,$om,$ol);
    my $da=$cosmo[5];
#    my $da=$dl/((1+$z)**2);
#    my $re_k=1454.44*$da*$re;
#    return ($ludis,$angsize,$propmo,$lbtime,$dm,$scale);
    my $re_k=$da*$re;
    return $re_k;
}


sub L_abs_lambda {
    my @data=@_;
    my $ho=$data[0];
    my $om=$data[1];
    my $ol=$data[2];
    my $z=$data[3];
    my $l=$data[4];
    my @cosmo=cosm($z,$ho,$om,$ol);

    my $dl=$cosmo[0];
    my $ratio=3.08567758e24; #Mpc to cm.
    $dl=$dl*$ratio;
    my $L=4*3.1416*($dl**2)*$l/(1+$z);
    return $L;
}


sub l_abs_lambda {
    my @data=@_;
    my $ho=$data[0];
    my $om=$data[1];
    my $ol=$data[2];
    my $z=$data[3];
    my $l=$data[4];
    my @cosmo=cosm($z,$ho,$om,$ol);

    my $dl=$cosmo[0];
#    print "$dl\n";
    my $ratio=3.08567758e24; #Mpc to cm.
    $dl=$dl*$ratio;
    my $L=4*3.1416*($dl**2)*$l/(1+$z);
    return $L;
}

sub QKS {
    local ($lambda)=@_;
    my $j;
    my $out=0;
    for ($j=1;$j<100;$j++) {
	$out=$out+2*((-1)**($j-1))*exp(-2*($j*$lambda)**2);
    }
    if ($out>1) {
	$out=1;
    }
    return $out;
}


sub KS_test_INT {
    my $n1=$_[0];
    my $d1=$_[1];
    my $n2=$_[2];
    my $d2=$_[3];
#    local ($n1,$d1,$n2,$d2)=@_;
    my @data1=@$d1;
    my @data2=@$d2;
    my @data_tmp;
    my @out;
    my $i;
    my $j;
    my @cdf1;
    my @cdf2;
    my $x;
    my $min1,$max1;
    my $min2,$max2;
#    my $a,$b;

    my $pdl_cdf1;
    my $pdl_cdf2;

    $n1=$#data1+1;
    $n2=$#data2+1;
    ($min1,$max1)=minmax(@data1);
    ($min2,$max2)=minmax(@data2);

    if ($min1>$min2) {
	$min1=$min2;
    }
    if ($max1<$max2) {
	$max1=$max2;
    }
    my @ord1 = sort {$a <=> $b} @data1;
    my @ord2 = sort {$a <=> $b} @data2;

    for ($i=0;$i<$n1;$i++) {
	if ($i==0) {
	    $cdf1[$i]=1/$n1;
	} else {
	    $cdf1[$i]=$cdf1[$i-1]+1/$n1;
	}	
    }
    for ($i=0;$i<$n2;$i++) {
	if ($i==0) {
	    $cdf2[$i]=1/$n2;
	} else {
	    $cdf2[$i]=$cdf2[$i-1]+1/$n2;
	}	
    }
    
    my $N=($n1+$n2);
    my @x=(0 .. $N);
    my $pdl_x = pdl(@x);
    $pdl_x = $min1 + (($max1-$min1)/$N)*$pdl_x;
#    print "\n VAL | $max1 | $min1 | $N | \n";
#    print "\n $pdl_x\n";

#    print "\n DATA1[$n1]= @data1 \n";
#    print "\n ORD1[$n1]= @ord1 \n";
#    print "\n CDF1[$n1]= @cdf1\n";

#    print "\n DATA2[$n2]= @data2 \n";
#    print "\n ORD2[$n2]= @ord2 \n";
#    print "\n CDF2[$n2]= @cdf2 \n";


#    print "PDL_x = $pdl_x\n";
    $pdl_cdf1 = interpol($pdl_x, pdl(@ord1) , pdl(@cdf1));
    $pdl_cdf2 = interpol($pdl_x, pdl(@ord2) , pdl(@cdf2));

    my $D=abs($pdl_cdf1-$pdl_cdf2);
    ($min2,$max2)=minmax(list($D));
    print "NNN = $N $n1 $n2 $min2 $max2\n";
    open(DIST,">ks_test.dist.txt");
    for ($i=0;$i<$N;$i++) {
	my $val1=$pdl_x->at($i);
	my $val2=$pdl_cdf1->at($i);
	my $val3=$pdl_cdf2->at($i);
	print DIST "$val1 $val2 $val3\n";
	print "$i $N $val1 $val2 $val3\n";
    }
    close(DIST);
    return $max2;
}


sub KS_test_BAD {
    local ($n1,$d1,$n2,$d2)=@_;
    my @data1=@$d1;
    my @data2=@$d2;
    my @data_tmp;
    my @out;
    my $i;
    my $j;
    #
    # First we order the data in @data1;
    #
    for ($j=0;$j<$n1;$j++) {
	$data_tmp[$j]=$data1[$j];
	for ($i=$j;$i<$n1;$i++) {
	    if ($data1[$i]<=$data_tmp[$j]) {
		$tmp=$data1[$i];
		$data1[$i]=$data_tmp[$j];
		$data_tmp[$j]=$tmp;
	    }
	}
    }
    #
    # Now  we look for D
    #
    my $d_max=0;
    my $d=0;
    my $cum1=0;
    my $cum2=0;
    for ($j=0;$j<$n1;$j++) {
	$cum1=$cum1+1/$n1;
	for ($i=0;$i<$n2;$i++) {
	    if ($data2[$i]<$data_tmp[$j]) {
		if ($j==0) {
		    $cum2=$cum2+1/$n2;
		} else {
		    if ($data2[$i]>=$data_tmp[$j-1]) {
			$cum2=$cum2+1/$n2;
		    }
		}
	    }
	}
	$d=abs($cum1-$cum2);
	if ($d_max<$d) {
	    $d_max=$d;
	}
    }
    return $d_max;
}




sub KS_test_2D {
    local ($n1,$dx1,$dy1,$n2,$dx2,$dy2)=@_;
    my @data_x1=@$dx1;
    my @data_y1=@$dy1;
    my @data_x2=@$dx2;
    my @data_y2=@$dy2;
    my @f1;
    my @f2;
    my $i;
    my $j;
    my $k;

    #
    # Now  we look for D
    #
    my $d_max=0;
    my $d=0;
    my $cum1=0;
    my $cum2=0;
    for ($j=0;$j<$n1;$j++) {
	$f2[0]=$f2[1]=$f2[2]=$f2[3]=0;
	$f1[0]=$f1[1]=$f1[2]=$f1[3]=0;
	$d=0;

	for ($i=0;$i<$n1;$i++) {
	    if (($data_x1[$i]<=$data_x1[$j])&&($data_y1[$i]>$data_y1[$j])) {
		$f1[0]=$f1[0]+1/$n1;
	    }
	    if (($data_x1[$i]>$data_x1[$j])&&($data_y1[$i]>$data_y1[$j])) {
		$f1[1]=$f1[1]+1/$n1;
	    }
	    if (($data_x1[$i]<=$data_x1[$j])&&($data_y1[$i]<=$data_y1[$j])) {
		$f1[2]=$f1[2]+1/$n1;
	    }
	    if (($data_x1[$i]>$data_x1[$j])&&($data_y1[$i]<=$data_y1[$j])) {
		$f1[3]=$f1[3]+1/$n1;
	    }
	}
	for ($i=0;$i<$n2;$i++) {
	    if (($data_x2[$i]<=$data_x1[$j])&&($data_y2[$i]>$data_y1[$j])) {
		$f2[0]=$f2[0]+1/$n2;
	    }
	    if (($data_x2[$i]>$data_x1[$j])&&($data_y2[$i]>$data_y1[$j])) {
		$f2[1]=$f2[1]+1/$n2;
	    }
	    if (($data_x2[$i]<=$data_x1[$j])&&($data_y2[$i]<=$data_y1[$j])) {
		$f2[2]=$f2[2]+1/$n2;
	    }
	    if (($data_x2[$i]>$data_x1[$j])&&($data_y2[$i]<=$data_y1[$j])) {
		$f2[3]=$f2[3]+1/$n2;
	    }
	   
	}
	for ($k=0;$k<4;$k++) {
	    if ($d<abs($f2[$k]-$f1[$k])) { 
		$d=abs($f2[$k]-$f1[$k]);
	    }
	}

	if ($d_max<$d) {
	    $d_max=$d;
	}	
    }
    return $d_max;
}



sub get_redshift_from_ned {

# Get parameters
    my $name = @_[0];

#    print "$name\n";

# contact the server
    if (open_TCP(F, "nedwww.ipac.caltech.edu", 80) == undef) {
	print "Error connecting to the NED server\n";
	exit(-1);
    }


    my $string="/cgi-bin/nph-datasearch?objname=".$name."&search_type=Redshifts&zv_breaker=30000.0&search_type=Redshifts";

#    print "\n$string\n";

    print F "GET $string HTTP/1.0\n";
    print F "Accept: */*\n";
    print F "User-Agent: hcat/1.0\n\n";


    my $hay_redshift=0;
    my $n_z=0;
    my $z=0;
    while ($linea=<F>) {

	#****************
	#print "$linea\n";

	if ($linea =~ "#No") {
	    $z[$n_z]=substr($linea,67,8);
#	    print "$n_z $z[$n_z]\n";
	    $n_z++;
	}
	
}
close(F);

if ($n_z==0) {
#    print "No hay redshift para el objeto $name\n";
    $z_tot=-1;
} else {

    if ($n_z==1) {
	$z_tot=$z[0];
#	print "El refshift es=$z[0]\n";
    } else {
#	print "Hay $n_z medidas del redshift:\n";
	$z_tot=0;
	for ($h=0;$h<$n_z;$h++) {
	    $z_tot=$z_tot+$z[$h];
#	    print "h=$h, z=$z[$h], $z_tot\n";
	}

	$z_tot=$z_tot/$n_z;
#	print "z_tot=$z_tot $n_z\n";
    }
}

    #print "$z_tot\n";
    return $z_tot;
}


sub get_data_sdss {
    my $r_ra = @_[0];
    my $r_dec = @_[1];
    my $url="http://cas.sdss.org/dr3/en/tools/search/x_radial.asp?ra=".$r_ra."+&dec=".$r_dec."&radius=0.02&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20&min_z=0&max_z=20&entries=top&topnum=10&format=csv";
    my $call="wget -q -r -l2 --no-parent -nd '$url' -O junk.html > /tmp/junk.txt ";
#    print "$url\n";
    system($call);
#  print "$call";
#  exit;

    open(WGET,"<junk.html");
    my  $out=<WGET>;
    my  $out=<WGET>;
    chop($out);
    my  @out_data=split(",",$out);	
    close(WGET);
    my  $call="rm junk.html";
    system($call);
    return @out_data;
}

sub get_data_ned {

# Get parameters
    my $name = @_[0];

    $name =~ s/\+/%2B/;
#    print "$name\n";
    $url="http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=".$name."&extend=no&img_stamp=YES";
#    $url="http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=".$name."&extend=no&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES";
    system("wget -q -r -l2 --no-parent -nd '$url' -O junk.html > /tmp/junk.txt ");
#   print "$url\n";
#   exit;
#    system("konqueror junk.html");
    my $hay_redshift=0;
    my $n_z=0;
    my $z=0;
    my $alpha=0;
    my $delta=0;
    my $ra=0;
    my $dec=0;
    open(F,"<junk.html");
    while ($linea=<F>) {
#	print "$linea\n";
	if (($linea =~ "Redshift")&&($linea =~ /\:/)) {
	    @data=split(" ",$linea);
	    $z=$data[2];
	}
	if (($linea =~ "Equatorial")&&($linea =~ "J2000.0")) {
	    @data=split(" ",$linea);
	    $alpha=$data[2];
	    $delta=$data[3];
	    $ra=$data[4];
	    $dec=$data[5];
	}
	
    }
    close(F);
    system("rm junk.html");

    return ($z,$alpha,$delta,$ra,$dec);
}



sub delay {
    my $delta=@_;
    $initial_time=time;
    do {
	$actual_time=time-$initial_time;
    } while ($delta<($actual_time));

}

sub open_TCP
{
    # get parameters
    my ($FS, $dest, $port) = @_;
    
    my $proto = getprotobyname('tcp');
    socket($FS, PF_INET, SOCK_STREAM, $proto);
    my $sin = sockaddr_in($port,inet_aton($dest));
    connect($FS,$sin) || return undef;
    
    my $old_fh = select($FS); 
    $| = 1;                       # don't buffer output
    select($old_fh);
    1;
}

sub mag_suma {
  local($mag1,$mag2)=@_;
  my $mag=-2.5*log10(10**(-0.4*$mag1)+10**(-0.4*$mag2));
  return $mag;
}

sub mag_resta {
  local($mag1,$mag2)=@_;
  my $mag=-2.5*log10(10**(-0.4*$mag1)-10**(-0.4*$mag2));
  return $mag;
}

sub ratio_mag {
  local($mag1,$mag2)=@_;
  my $ratio=10**(-0.4*($mag2-$mag1));
  return $ratio;
}


sub gauss_filter {
    my $sigma=$_[0];
    my $box=3*int($sigma);
    my $array=$_[1];
    my $norm=($sigma*sqrt(2*3.1416));
  
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    my @val=@$array;
    my $i,$jk;
    for ($i=$box;$i<($#val+1-$box);$i++) {
	$val[$i]=0;
	for ($jk=0;$jk<2*$box;$jk++) {
	    $tmp=@$array[$i-$box+$jk]*(exp(-0.5*(($jk-$box)/$sigma)**2))/$norm;
	    $val[$i]=$val[$i]+$tmp;
	}
#	$val[$i]=mean(@tmp);
    }

    for ($i=1;$i<$box;$i++) {
	my @tmp;
	$effec_box=$i;
	for ($jk=0;$jk<2*$effec_box;$jk++) {
	    $tmp[$jk]=@$array[$i-$effec_box+$jk];
	}
	$val[$i]=mean(@tmp);
#	$val[$i]=$val[$box+1];
    }

    for ($i=($#val+1-$box);$i<$#val;$i++) {
	my @tmp;
	$effec_box=$#val-$i+1;
	for ($jk=0;$jk<2*$effec_box;$jk++) {
	    $tmp[$jk]=@$array[$i-$effec_box+$jk];
	}
	$val[$i]=mean(@tmp);
#	$val[$i]=$val[$#val-$box-1];
    }
    return @val;
}

sub pgpoint3d {
    my $xmin=$_[0];
    my $xmax=$_[1];
    my $ymin=$_[2];
    my $ymax=$_[3];
    my $zmin=$_[4];
    my $zmax=$_[5];
    my $angle=($_[6])*3.1416/180;
    my @x=@$_[7];
    my @y=@$_[8];
    my @z=@$_[9];
    pgsvp(0.17,0.83,0.15,0.95);
    pgswin(0,1+sin($angle),0,1+cos($angle));
    my $i;
    for ($i=0;$i<=$#x;$i++) {
	$x[$i]=($x[$i]-$xmin)/($xmax-$xmin);
	$y[$i]=($y[$i]-$ymin)/($ymax-$ymin);
	$z[$i]=($z[$i]-$ymin)/($ymax-$ymin);
    }

    return;
}

#
# A_l/A_v as a function of lambda and Rv
# obtained from Cardielli, Clayton & Mathis, 1989, ApJ, 345, 245
#
sub A_l {
    my $Rv=$_[0];    
    my $l=$_[1];
    my $l=$l/10000; #Amstrongs to Microns
#    if ($l<1) {
    my $x;
    my $y;
    my $ax;
    my $bx;
 
    $x=1/$l;
    if ($x>1.1) {
	$y=($x-1.82);
	$ax=1+0.17699*$y-0.50447*$y**2-0.02427*$y**3+0.72085*$y**4+0.01979*$y**5-0.77530*$y**6+0.32999*$y**7;
	$bx=1.41338*$y+2.28305*$y**2+1.07233*$y**3-5.38434*$y**4-0.62251*$y**5+5.30260*$y**6-2.09002*$y**7;
    } else {
	$ax=0.574*$x**1.61;
	$bx=-0.527*$x**1.61;
    }

	$Arat=$ax+$bx/$Rv;
#    } else {
#	my $Arat=1;
#    }
#    if ($l>=1.2) {
#	$Arat=0.26*$x;
#    }
#    print "$l ";
    return $Arat;
}

sub A_l3 {
    my $Rv=$_[0];    
    my $l=$_[1];
    my $l=$l/10000; #Amstrongs to Microns
#    if ($l<1) {
    my $x;
    my $y;
    my $ax;
    my $bx;
 
    $x=1/$l;
    
    $Arat=1+(0.404/$Rv)*$x;
    if ($x>1.1) {
	$y=($x-1.82);
	$ax=1+0.17699*$y-0.50447*$y**2-0.02427*$y**3+0.72085*$y**4+0.01979*$y**5-0.77530*$y**6+0.32999*$y**7;
	$bx=1.41338*$y+2.28305*$y**2+1.07233*$y**3-5.38434*$y**4-0.62251*$y**5+5.30260*$y**6-2.09002*$y**7;
    } else {
	$ax=0.574*$x**1.61;
	$bx=-0.527*$x**1.61;
    }

#	$Arat=$ax+$bx/$Rv;
#    } else {
#	my $Arat=1;
#    }
#    if ($l>=1.2) {
#	$Arat=0.26*$x;
#    }
#    print "$l ";
    return $Arat;
}


sub A_l2 {
    # Calzetti
    my $Rv=$_[0];    
    my $l=$_[1];
    $l=$l/10000; #Amstrongs to Microns
    my $x=1/$l;
    my $kl=-2.156+1.509*$x-0.198*$x**2+0.011*$x**3;
    if ($kl>0) {
	my $Arat=$kl/$Rv;
    } else {
	my $Arat=0;
    }
    return $Arat;
}



1;

sub probks {
    my $alam=$_[0];
    my $j;
    my $a2;
    my $fac=2;
    my $sum=0;
    my $term;
    my $termbf=0;
    my $EPS1=0.001;
    my $EPS2=1e-8;
    $a2=-2*$alam*$alam;
    for ($j=0;$j<100;$j++) {
	$term=$fac*exp($a2*$j*$j);
	$sum=$sum+$term;
	if ((abs($term)<= ($EPS1*$termbf))||(abs($term)<= ($EPS2*$sum))) {
	    return $sum;
	}
	$fac=(-1)*$fac;
	$termbf=abs($term);
    }
    return 1;
}



sub KS_test {
    my $n1=$_[0];
    my $d1=$_[1];
    my $n2=$_[2];
    my $d2=$_[3];
#    local ($n1,$d1,$n2,$d2)=@_;
    my @data1=@$d1;
    my @data2=@$d2;
    my @data_tmp;
    my @out;
    my $i;
    my $j;
    my @cdf1;
    my @cdf2;
    my $x;
    my $min1,$max1;
    my $min2,$max2;
#    my $a,$b;

    my $pdl_cdf1;
    my $pdl_cdf2;

    $n1=$#data1+1;
    $n2=$#data2+1;
    ($min1,$max1)=minmax(@data1);
    ($min2,$max2)=minmax(@data2);

    if ($min1>$min2) {
	$min1=$min2;
    }
    if ($max1<$max2) {
	$max1=$max2;
    }
    my @ord1 = sort {$a <=> $b} @data1;
    my @ord2 = sort {$a <=> $b} @data2;

    for ($i=0;$i<$n1;$i++) {
	if ($i==0) {
	    $cdf1[$i]=1/$n1;
	} else {
	    $cdf1[$i]=$cdf1[$i-1]+1/$n1;
	}	
    }
    for ($i=0;$i<$n2;$i++) {
	if ($i==0) {
	    $cdf2[$i]=1/$n2;
	} else {
	    $cdf2[$i]=$cdf2[$i-1]+1/$n2;
	}	
    }
    
    my $N=2*($n1+$n2);
    my @x=(0 .. $N);
    my $pdl_x = pdl(@x);
    $pdl_x = $min1 + (($max1-$min1)/$N)*$pdl_x;
    $pdl_cdf1=zeroes($N);
    $pdl_cdf2=zeroes($N);
    for ($i=0;$i<$n1;$i++) {
	my $val=$ord1[$i];
#	print "MY $J=int(($val-$min1)/(($max1-$min1)/$N))-1\n";
	my $J;
	if ($max1!=$min1) {
	    $J=int(($val-$min1)/(($max1-$min1)/$N))-1;
	} else {
	    $J=$N;
	}
	for ($j=$J;$j<$N;$j++) {
	    my $val=$pdl_cdf1->at($j);
	    $val=$val+1/$n1;
	    set($pdl_cdf1,$j,$val);
	}
    }
    for ($i=0;$i<$n2;$i++) {
	my $val=$ord2[$i];
	my $J;
	if ($max1!=$min1) {
	    $J=int(($val-$min1)/(($max1-$min1)/$N))-1;
	} else {
	    $J=$N;
	}
	for ($j=$J;$j<$N;$j++) {
	    my $val=$pdl_cdf2->at($j);
	    $val=$val+1/$n2;
	    set($pdl_cdf2,$j,$val);
	}
    }

#    print "$pdl_cdf1\n $pdl_cdf2\n";
#    $pdl_cdf1 = interpol($pdl_x, pdl(@ord1) , pdl(@cdf1));
#    $pdl_cdf2 = interpol($pdl_x, pdl(@ord2) , pdl(@cdf2));

    my $D=abs($pdl_cdf1-$pdl_cdf2);
    ($min2,$max2)=minmax(list($D));
    open(DIST,">ks_test.dist.txt");
    for ($i=0;$i<$N;$i++) {
	my $val1=$pdl_x->at($i);
	my $val2=$pdl_cdf1->at($i);
	my $val3=$pdl_cdf2->at($i);
	print DIST "$val1 $val2 $val3\n";
    }
    close(DIST);
    return $max2;
}


sub cal_ctab() {
    $nc=254;
    open(CTAB,"</home/sanchez/sda2/code/R3D/cal_ctab.txt");
    while($cmap=<CTAB>) {
	chop($cmap);
	@data=split(" ",$cmap);
	$nc=$data[0];
#    $nc++;
	$r[$nc-1]=$data[1]/255;
	$g[$nc-1]=$data[2]/255;
	$b[$nc-1]=$data[3]/255;
	$l[$nc]=$nc/255;
#		print "";
#    $cmap=<DATA>
    }
    close(CTAB);
    $bright=1; $contrast=0.5;
    $r[0]=1.0;
    $g[0]=1.0;
    $b[0]=1.0;
    $nc=254;
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);
      
}

sub cal_ctab_bc() {
    $nc=254;
    open(CTAB,"</home/sanchez/sda2/code/R3D/cal_ctab.txt");
    while($cmap=<CTAB>) {
	chop($cmap);
	@data=split(" ",$cmap);
	$nc=$data[0];
#    $nc++;
	$r[$nc-1]=$data[1]/255;
	$g[$nc-1]=$data[2]/255;
	$b[$nc-1]=$data[3]/255;
	$l[$nc]=$nc/255;
#		print "";
#    $cmap=<DATA>
    }
    close(CTAB);
#    $bright=1; $contrast=0.5;
     if ($bright>0) {
	 $r[0]=1.0;
	 $g[0]=1.0;
	 $b[0]=1.0;
     } else {
	 $r[253]=1.0;
	 $g[253]=1.0;
	 $b[253]=1.0;
     }



    $nc=254;
    pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);
      
}

sub cal_vel_ctab() {
     open(CTAB,"</home/sanchez/sda2/code/R3D/cal_vel_ctab.txt");
     $nc0=255;
     while($cmap=<CTAB>) {
	 chop($cmap);
	 @data=split(" ",$cmap);
	 $nc=int($data[0]);
	 $k=256-$nc-1;
	 $r[$k]=$data[1]/255;
	 $g[$k]=$data[2]/255;
	 $b[$k]=$data[3]/255;
	 
	 $k=$nc-1;
	 $rr[$k]=$data[1]/255;
	 $gg[$k]=$data[2]/255;
	 $bb[$k]=$data[3]/255;
	 
	 $l[$nc-1]=$nc/255;
     }
     close(CTAB);
     $bright=1; $contrast=0.5;
     pgscir(50,128+50); $nc=254;
#pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);
	 $rr[0]=1.0;
	 $gg[0]=1.0;
	 $bb[0]=1.0;
	 $r[0]=1.0;
	 $g[0]=1.0;
	 $b[0]=1.0;    
     pgctab(\@l,\@rr,\@gg,\@bb,$nc,$bright,$contrast);
}


sub cal_vel_ctab_bc() {
     open(CTAB,"</home/sanchez/sda2/code/R3D/cal_vel_ctab.txt");
     $nc0=255;
     while($cmap=<CTAB>) {
	 chop($cmap);
	 @data=split(" ",$cmap);
	 $nc=int($data[0]);
	 $k=256-$nc-1;
	 $r[$k]=$data[1]/255;
	 $g[$k]=$data[2]/255;
	 $b[$k]=$data[3]/255;
	 
	 $k=$nc-1;
	 $rr[$k]=$data[1]/255;
	 $gg[$k]=$data[2]/255;
	 $bb[$k]=$data[3]/255;
	 
	 $l[$nc-1]=$nc/255;
     }
     close(CTAB);
#     $bright=1; $contrast=0.5;
     pgscir(50,128+50); $nc=254;
#pgctab(\@l,\@r,\@g,\@b,$nc,$bright,$contrast);
     if ($bright>0) {
	 $rr[0]=1.0;
	 $gg[0]=1.0;
	 $bb[0]=1.0;
	 $r[0]=1.0;
	 $g[0]=1.0;
	 $b[0]=1.0;    
} else {
	 $rr[253]=1.0;
	 $gg[253]=1.0;
	 $bb[253]=1.0;
	 $r[253]=1.0;
	 $g[253]=1.0;
	 $b[253]=1.0;    
} 
     pgctab(\@l,\@rr,\@gg,\@bb,$nc,$bright,$contrast);
}
