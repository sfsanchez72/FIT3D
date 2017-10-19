#########################################
# My Subroutines
# S.F.Sanchez 2004
##########################################


#$naxes=read_naxes($fitsfile);
#@in_data=read_img($fitsfile);
#$npix=@$naxes[0];
#$nb_spec=@$naxes[1];
#write_img($outfile,$npix,$nb_spec,\@out_data); 

use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Primitive;

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
    for ($i=0;$i<$#images;$i++) {
	my @aa_tmp=read_img($images[$i]);
	for ($j=0;$j<$nb_spec;$j++) {
	    for ($k=0;$k<$npix;$k++) {
		for ($i=0;$i<$#images;$i++) {
		    $a_tmp[$j][$k][$i]=$aa_tmp[$j][$k];
		}
	    }
	}
    }

    for ($j=0;$j<$nb_spec;$j++) {
	for ($k=0;$k<$npix;$k++) {
	    my @f_tmp;
	    for ($i=0;$i<$#images;$i++) {
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
    my $n=$#wa;
    my $E_BV=$_[1];
#    my @E_L=(1000,1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000);
#    my @E_LV=(10.0,9.1,8.3,7.7,7.1,6.7,6.2,5.6,5.0,4.5,4.2,3.8,3.6,3.3,2.9,2.5);
    my $x=pdl(40000,26500,12200,6000,5470,4670,4110,2700,2600);
    my $y=pdl(0.000,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591);
    my $yi = interpol(pdl($w), $x, $y);
    my $yi = $yi;
    my $j,@out;
    for ($j=0;$j<$n;$j++) {
	$out[$j] = $E_BV*$yi->slice($j)->sclr;
    }
    return @out;
}

sub median {
  local(@data)=@_;
# sort numerically ascending
  my @out_data = sort {$a <=> $b} @data;
  my $n_median=int($#out_data*0.5);
  my $median = ($out_data[$n_median-1]+$out_data[$n_median]+$out_data[$n_median+1])/3;
  return $median;
}

sub median_filter {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    @val=@$array;
    my $i,$j;
    for ($i=$box;$i<($#val-$box);$i++) {
	my @tmp;
	for ($j=0;$j<2*$box;$j++) {
	    $tmp[$j]=@$array[$i-$box+$j];
	}
	$val[$i]=median(@tmp);
    }

    for ($i=1;$i<$box;$i++) {
	my @tmp;
	$effec_box=$i;
	for ($j=0;$j<2*$effec_box;$j++) {
	    $tmp[$j]=@$array[$i-$effec_box+$j];
	}
	$val[$i]=median(@tmp);
#	$val[$i]=$val[$box+1];
    }

    for ($i=($#val-$box);$i<($#val-1);$i++) {
	my @tmp;
	$effec_box=$#val-$i;
	for ($j=0;$j<2*$effec_box;$j++) {
	    $tmp[$j]=@$array[$i-$effec_box+$j];
	}
	$val[$i]=median(@tmp);
#	$val[$i]=$val[$#val-$box-1];
    }
    return @val;
}


sub median_box {
    my $box=$_[0];
    my $array=$_[1];
    if ($box=(2*int($box/2))) {
	$box=$box+1;
    }
    @in_val=@$array;
    my $i,$j,$k;
    $k=0;
    for ($i=$box;$i<($#in_val-$box);$i=$i+2*$box) {
	my @tmp;
	for ($j=0;$j<2*$box;$j++) {
	    $tmp[$j]=@$array[$i-$box+$j];
	}
	$out_val[$k]=median(@tmp);
#	$out_sval[$k]=sigma(@tmp);
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
    for ($i=$box;$i<($#in_val-$box);$i=$i+2*$box) {
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


sub mean {
  local(@data)=@_;
  my $sum=0; 
  my $j;
  for ($j=0;$j<$#data;$j++) {
      $sum=$sum+$data[$j];
  }
  if ($#data>0) {
      my $mean = $sum/$#data;
  }
  return $mean;
}

sub sum_over {
  local(@data)=@_;
  my $sum=0; 
  my $j;
  for ($j=0;$j<$#data;$j++) {
      $sum=$sum+$data[$j];
  }
  return $sum;
}

sub sigma {
  local(@data)=@_;
  my $mean = mean(@data);
  my $stddev = 0;
  my $j;
  my $sum=0; 
  for ($j=0;$j<$#data;$j++) {
      $sum=$sum+($data[$j]-$mean)**2;
  }
  if ($#data>1) {
      $sum=$sum/($#data-1);
  }
  $stddev=sqrt($sum);
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

sub read_img_headers() {
    my $file=$_[0];
    my $tmp=$_[1];
    my @headers_id=@$tmp;
    my $status = 0;
    my $fptr = Astro::FITS::CFITSIO::open_file($file,Astro::FITS::CFITSIO::READONLY(),$status);
    my $naxes;
    my $value,$comment;
    my @values;
    my $i;
    for ($i=0;$i<$#headers_id+1;$i++) {
	my $header_id=$headers_id[$i];
	$fptr->read_keyword($header_id,$value,$comment,$status);
	$values[$i]=$value;
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
     $out->write_key_flt('CRVAL1',$start,5,'CRVAL',$status) and print "CRVAL status = $status\n";
     $out->write_key_lng('CRPIX1',1,'CRpix1',$status) and print "CRPIX1 status = $status\n";
     $out->write_key_flt('CD1_1',$delta,5,'Cd1_1',$status) and print "CD1_1 status = $status\n";
     $out->write_key_flt('CDELT1',$delta,5,'Cdelta1',$status) and print "CDELT1 status = $status\n";
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




sub log10 {
  my $n = shift;
  if ($n==0) {
      $n=1;
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
	for ($j=0;$j<$#data;$j++) {
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

    my $ludis=$r*((1+$z)**2);
    my $angsize=$r;
    my $propmo=$r*(1+$z);
    my $lbtime=$th*tint($z,$om,$ol);
    my $dm=5*0.4343*log($r*((1+z)**2))+25;
    my $scale=$angsize/206.264806;

    $cosmo->{ludis}=$ludis;
    $cosmo->{angsize}=$angsize;
    $cosmo->{propmo}=$propmo;
    $cosmo->{lbtime}=$lbtime;
    $cosmo->{dm}=$dm;
    $cosmo->{scale}=$scale;
    
    return ($ludis,$angsize,$propmo,$lbtime,$dm,$scale);
}



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




1;

