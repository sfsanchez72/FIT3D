sub shift_convolve {
    my $junk=$_[0];
    my @wave=@$junk;
    my $junk=$_[1];
    my @flux=@$junk;
    my $redshift=$_[2];
    my $sigma_inst=$_[3];
    my $sigma_km_h=$_[4];
    my $dpix=$wave[1]-$wave[0];
    my $rsigma=$sigma_inst/$dpix;
    my $pdl_flux=pdl(@flux);
    my $e=exp(1);
    my $pdl_wave=pdl(@wave);
######################################################
# Convolved with Sigma_inst
#
    my $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    my $kernel=zeroes(2*$box+1);
    my    $norm=0;
    for ($j=0;$j<2*$box+1;$j++) {
	my $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
}
    $kernel=$kernel/$norm;
    my $pdl_flux_conv_inst = conv2d $pdl_flux,$kernel;
#
#
###################################################
    my $l_pdl_wave=log($pdl_wave);
    my $l_pdl_wave0=$l_pdl_wave->at(0);
    my $l_pdl_wave1=$l_pdl_wave->at($#wave);
# Take care: Factor 5!!!!
# Interpolation Factor
my $f_fine=3;
    my $d_l_pdl_wave=($l_pdl_wave->at(1)-$l_pdl_wave->at(0))/$f_fine;
    my $n_l_pdl_wave=int(($l_pdl_wave1-$l_pdl_wave0)/($d_l_pdl_wave));
    my $pdl_ones=pdl(0..$n_l_pdl_wave);
    my $l_pdl_wave_l=$l_pdl_wave0+$d_l_pdl_wave*$pdl_ones;
    my $pdl_wave_l=$e**($l_pdl_wave_l);
    my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);

######################################################
# Convolved with Sigma_km_h
#
    my $rsigma_km_h=($sigma_km_h/300000)/($d_l_pdl_wave);#*$f_fine);
    my $box=int(3*$rsigma_km_h);
    if ($box<3) {
	$box=3;
    }
    my $kernel=zeroes(2*$box+1);
    my    $norm=0;
for ($j=0;$j<2*$box+1;$j++) {
    my $gaus=exp(-0.5*((($j-$box)/$rsigma_km_h)**2));    
    set($kernel,$j,$gaus);
    $norm=$norm+$gaus;
}
    $kernel=$kernel/$norm;
    my $pdl_flux_conv_km_h = conv2d $pdl_flux_conv_inst,$kernel;
#
#
###################################################
    my $pdl_wave_l_redshift=$pdl_wave_l/(1+$redshift);
    my $pdl_wave_redshift=$pdl_wave*(1+$redshift);
#    my $pdl_flux_l_wave_km_h=interpol($pdl_wave_l_redshift,$pdl_wave,$pdl_flux_conv_km_h);

    my $pdl_flux_l_wave_km_h=interpol($pdl_wave,$pdl_wave_redshift,$pdl_flux_conv_km_h);
    return $pdl_flux_l_wave_km_h;

}


sub shift_convolve_pdl_WW1 {
#     ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
#     print "start shift_convolve_pdl\n";
    my $pdl_wave_out=$_[0];
    my $pdl_wave_in=$_[1];
    my $pdl_flux_in=$_[2];
    my $redshift=$_[3];
    my $sigma_inst=$_[4];
    my $sigma_km_h=$_[5];
    $sigma_inst=$sigma_inst*(1+$redshift);
    $sigma_km_h=$sigma_km_h/(1+$redshift);
    my $sub_pix=1+rand();
    my ($wmin,$wmax)=minmax(list($pdl_wave_in));
    my $dpix_ini=$pdl_wave_in->at(1)-$pdl_wave_in->at(0);
    my $dpix_out=$pdl_wave_out->at(1)-$pdl_wave_out->at(0);

    my $dpix=$dpix_ini/$sub_pix;
    my $n_sub=($wmax-$wmin)/$dpix;
    my $pdl_wave=$wmin+(($wmax-$wmin)/$n_sub)*pdl(0..$n_sub);
    my $pdl_flux=interpol($pdl_wave,$pdl_wave_in,$pdl_flux_in);


    my $rsigma=$sigma_inst/$dpix;
    my $e=exp(1);

#    my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
######################################################
# Convolved with Sigma_inst
#
#    my $box=int(3*$rsigma);

    my $pdl_flux_conv_inst;
    if ($sigma_inst>0.1) {
	my $box=int(3*$sigma_inst*2/$dpix);
	if ($box<3) {
	    $box=3;
	}
	my $kernel=zeroes(2*$box+1);
	my    $norm=0;
	for ($j=0;$j<2*$box+1;$j++) {
	    my $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	    set($kernel,$j,$gaus);
	    $norm=$norm+$gaus;
	}
	$kernel=$kernel/$norm;
	$pdl_flux_conv_inst = conv2d $pdl_flux,$kernel;
    } else {
	$pdl_flux_conv_inst=$pdl_flux;
    }

     
#
#
###################################################
    my $l_pdl_wave=log($pdl_wave);
    my $l_pdl_wave0=$l_pdl_wave->at(0);
    my $l_pdl_wave1=$l_pdl_wave->at($#wave);
# Take care: Factor 5!!!!
# Interpolation Factor
#my $f_fine=5;
#my $f_fine=5;
     my $f_fine=5;
    my $d_l_pdl_wave=($l_pdl_wave->at(1)-$l_pdl_wave->at(0))/$f_fine;
    my $n_l_pdl_wave=int(($l_pdl_wave1-$l_pdl_wave0)/($d_l_pdl_wave));
    my $pdl_ones=pdl(0..$n_l_pdl_wave);
    my $l_pdl_wave_l=$l_pdl_wave0+$d_l_pdl_wave*$pdl_ones;
    my $pdl_wave_l=$e**($l_pdl_wave_l);
    my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
######################################################
# Convolved with Sigma_km_h
#
    if ($sigma_km_h==0) {
	 $sigma_km_h=1;
     }
 
#   my $rsigma_km_h=($sigma_km_h/300000)/($d_l_pdl_wave);#/$dpix;
   my $rsigma_km_h=(($sigma_km_h/300000)/($d_l_pdl_wave));#/$dpix_out;
#    print "   my $rsigma_km_h=(($sigma_km_h/300000)/($d_l_pdl_wave))*$dpix;\n";

 #     print "    my $rsigma_km_h=($sigma_km_h/300000)/($d_l_pdl_wave*$f_fine)\n";
#    my $box=int(3*$rsigma_km_h);
    my $box=int(5*(500/300000)/($d_l_pdl_wave));
#     print "BOX,RSIGMA = $box $rsigma_km_h $d_l_pdl_wave\n";
    if ($box<3) {
	$box=3;
    }
    my $kernel=zeroes(2*$box+1);
    my    $norm=0;
#     print "#################\n";
    for ($j=0;$j<2*$box+1;$j++) {
	my $gaus=exp(-0.5*((($j-$box)/$rsigma_km_h)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
#	print "$j $gaus\n";
   }
#     print "#################\n";
    $kernel=$kernel/$norm;
#
# Check THIS!!! not wave_l???
#     my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
#    my $pdl_flux_conv_km_h = conv2d $pdl_flux_conv_inst,$kernel;
    my $pdl_flux_conv_km_h = conv2d $pdl_flux_l_wave,$kernel;
#
#
###################################################
    my $pdl_wave_l_redshift=$pdl_wave_l/(1+$redshift);
    my $pdl_wave_redshift=$pdl_wave*(1+$redshift);
    my $pdl_wave_out_redshift=$pdl_wave_out*(1+$redshift);
#    my $pdl_flux_l_wave_km_h=interpol($pdl_wave_l_redshift,$pdl_wave,$pdl_flux_conv_km_h);
 #   print "PASO\n";

#          print "Before interp\n";

    my $pdl_flux_l_wave_km_h=interpol($pdl_wave_out,$pdl_wave_l_redshift,$pdl_flux_conv_km_h);
#     print "After interp\n";
#     ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
#     print "end shift_convolve_pdl\n";
#
    return $pdl_flux_l_wave_km_h;
# 	 return $pdl_flux_conv_inst;   
}

sub shift_convolve_pdl {
#     ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
#     print "start shift_convolve_pdl\n";
    my $pdl_wave_out=$_[0];
    my $pdl_wave_in=$_[1];
    my $pdl_flux_in=$_[2];
    my $redshift=$_[3];
    my $sigma_inst=$_[4];
    my $sigma_km_h=$_[5];
    $sigma_inst=$sigma_inst*(1+$redshift);
    $sigma_km_h=$sigma_km_h/(1+$redshift);
    my $sub_pix=1+rand();
    my ($wmin,$wmax)=minmax(list($pdl_wave_in));
    my $dpix_ini=$pdl_wave_in->at(1)-$pdl_wave_in->at(0);
    my $dpix_out=$pdl_wave_out->at(1)-$pdl_wave_out->at(0);

    my $dpix=$dpix_ini/$sub_pix;
    my $n_sub=($wmax-$wmin)/$dpix;
    my $pdl_wave=$wmin+(($wmax-$wmin)/$n_sub)*pdl(0..$n_sub);
    my $pdl_flux=interpol($pdl_wave,$pdl_wave_in,$pdl_flux_in);


    my $rsigma=$sigma_inst/$dpix;
    my $e=exp(1);

#    my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
######################################################
# Convolved with Sigma_inst
#
#    my $box=int(3*$rsigma);

    my $pdl_flux_conv_inst;
    if ($sigma_inst>0.1) {
	my $box=int(3*$sigma_inst*2/$dpix);
	if ($box<3) {
	    $box=3;
	}
	my $kernel=zeroes(2*$box+1);
	my    $norm=0;
	for ($j=0;$j<2*$box+1;$j++) {
	    my $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	    set($kernel,$j,$gaus);
	    $norm=$norm+$gaus;
	}
	$kernel=$kernel/$norm;
	$pdl_flux_conv_inst = conv2d $pdl_flux,$kernel;
    } else {
	$pdl_flux_conv_inst=$pdl_flux;
    }

     
#
#
###################################################

    my $l_pdl_wave=log($pdl_wave);
    my $l_pdl_wave0=$l_pdl_wave->at(0);
    my $l_pdl_wave1=$l_pdl_wave->at($#wave);
# Take care: Factor 5!!!!
# Interpolation Factor
#my $f_fine=5;
#my $f_fine=5;
     my $f_fine=5;
    my $d_l_pdl_wave=($l_pdl_wave->at(1)-$l_pdl_wave->at(0))/$f_fine;
    my $n_l_pdl_wave=int(($l_pdl_wave1-$l_pdl_wave0)/($d_l_pdl_wave));
    my $pdl_ones=pdl(0..$n_l_pdl_wave);
    my $l_pdl_wave_l=$l_pdl_wave0+$d_l_pdl_wave*$pdl_ones;
    my $pdl_wave_l=$e**($l_pdl_wave_l)/(1+$redshift);
    my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
######################################################
# Convolved with Sigma_km_h
#
    if ($sigma_km_h==0) {
	 $sigma_km_h=1;
     }
 
#   my $rsigma_km_h=($sigma_km_h/300000)/($d_l_pdl_wave);#/$dpix;
   my $rsigma_km_h=(($sigma_km_h/300000)/($d_l_pdl_wave));#/$dpix_out;
#    print "   my $rsigma_km_h=(($sigma_km_h/300000)/($d_l_pdl_wave))*$dpix;\n";

 #     print "    my $rsigma_km_h=($sigma_km_h/300000)/($d_l_pdl_wave*$f_fine)\n";
#    my $box=int(3*$rsigma_km_h);
    my $box=int(5*(500/300000)/($d_l_pdl_wave));
#     print "BOX,RSIGMA = $box $rsigma_km_h $d_l_pdl_wave\n";
    if ($box<3) {
	$box=3;
    }
    my $kernel=zeroes(2*$box+1);
    my    $norm=0;
#     print "#################\n";
    for ($j=0;$j<2*$box+1;$j++) {
	my $gaus=exp(-0.5*((($j-$box)/$rsigma_km_h)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
#	print "$j $gaus\n";
   }
#     print "#################\n";
    $kernel=$kernel/$norm;
#
# Check THIS!!! not wave_l???
#     my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
#    my $pdl_flux_conv_km_h = conv2d $pdl_flux_conv_inst,$kernel;
    my $pdl_flux_conv_km_h = conv2d $pdl_flux_l_wave,$kernel;
#
#
###################################################
    my $pdl_wave_l_redshift=$pdl_wave_l;#/(1+$redshift);
#    my $pdl_wave_redshift=$pdl_wave*(1+$redshift);
#    my $pdl_wave_out_redshift=$pdl_wave_out*(1+$redshift);
#    my $pdl_flux_l_wave_km_h=interpol($pdl_wave_l_redshift,$pdl_wave,$pdl_flux_conv_km_h);
 #   print "PASO\n";

#          print "Before interp\n";

    my $pdl_flux_l_wave_km_h=interpol($pdl_wave_out,$pdl_wave_l_redshift,$pdl_flux_conv_km_h);
#     print "After interp\n";
#     ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
#     print "end shift_convolve_pdl\n";
#






    return $pdl_flux_l_wave_km_h;
# 	 return $pdl_flux_conv_inst;   
}

sub shift_convolve_pdl_NO {
#     ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
#     print "start shift_convolve_pdl\n";
    my $pdl_wave_out=$_[0];
    my $pdl_wave_in=$_[1];
    my $pdl_flux_in=$_[2];
    my $redshift=$_[3];
    my $sigma_inst=$_[4];
    my $sigma_km_h=$_[5];
    my $sub_pix=1+rand();
    my ($wmin,$wmax)=minmax(list($pdl_wave_in));
    my $dpix_ini=$pdl_wave_in->at(1)-$pdl_wave_in->at(0);

    my $dpix=$dpix_ini/$sub_pix;
    my $n_sub=($wmax-$wmin)/$dpix;
    my $pdl_wave=$wmin+(($wmax-$wmin)/$n_sub)*pdl(0..$n_sub);
    my $pdl_flux=interpol($pdl_wave,$pdl_wave_in,$pdl_flux_in);


    my $rsigma=$sigma_inst/$dpix;
    my $e=exp(1);

#    my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);


    my $pdl_flux_conv_inst=$pdl_flux;


    my $l_pdl_wave=log($pdl_wave);
    my $l_pdl_wave0=$l_pdl_wave->at(0);
    my $l_pdl_wave1=$l_pdl_wave->at($#wave);
# Take care: Factor 5!!!!
# Interpolation Factor
#my $f_fine=5;
#my $f_fine=5;
     my $f_fine=5;
    my $d_l_pdl_wave=($l_pdl_wave->at(1)-$l_pdl_wave->at(0))/$f_fine;
    my $n_l_pdl_wave=int(($l_pdl_wave1-$l_pdl_wave0)/($d_l_pdl_wave));
    my $pdl_ones=pdl(0..$n_l_pdl_wave);
    my $l_pdl_wave_l=$l_pdl_wave0+$d_l_pdl_wave*$pdl_ones;
    my $pdl_wave_l=$e**($l_pdl_wave_l);
    my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
######################################################
# Convolved with Sigma_km_h
#
    if ($sigma_km_h==0) {
	 $sigma_km_h=1;
     }
 
#   my $rsigma_km_h=($sigma_km_h/300000)/($d_l_pdl_wave);#/$dpix;
   my $rsigma_km_h=(($sigma_km_h/300000)/($d_l_pdl_wave));#*$dpix;
 #     print "    my $rsigma_km_h=($sigma_km_h/300000)/($d_l_pdl_wave*$f_fine)\n";
#    my $box=int(3*$rsigma_km_h);
    my $box=int(5*(500/300000)/($d_l_pdl_wave));
#     print "BOX,RSIGMA = $box $rsigma_km_h $d_l_pdl_wave\n";
    if ($box<3) {
	$box=3;
    }
    my $kernel=zeroes(2*$box+1);
    my    $norm=0;
#     print "#################\n";
    for ($j=0;$j<2*$box+1;$j++) {
	my $gaus=exp(-0.5*((($j-$box)/$rsigma_km_h)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
#	print "$j $gaus\n";
   }
#     print "#################\n";
    $kernel=$kernel/$norm;
#
# Check THIS!!! not wave_l???
#     my $pdl_flux_l_wave=interpol($pdl_wave_l,$pdl_wave,$pdl_flux_conv_inst);
#    my $pdl_flux_conv_km_h = conv2d $pdl_flux_conv_inst,$kernel;
    my $pdl_flux_conv_km_h = conv2d $pdl_flux_l_wave,$kernel;
#
#
###################################################
    my $pdl_wave_l_redshift=$pdl_wave_l/(1+$redshift);
    my $pdl_wave_redshift=$pdl_wave*(1+$redshift);
    my $pdl_wave_out_redshift=$pdl_wave_out*(1+$redshift);
#    my $pdl_flux_l_wave_km_h=interpol($pdl_wave_l_redshift,$pdl_wave,$pdl_flux_conv_km_h);
 #   print "PASO\n";

#          print "Before interp\n";
    my $pdl_flux_l_wave_km_h=interpol($pdl_wave_out,$pdl_wave_l_redshift,$pdl_flux_conv_km_h);

######################################################
# Convolved with Sigma_inst
#
#    my $box=int(3*$rsigma);
    my $pdl_flux=$pdl_flux_l_wave_km_h;
     my $pdl_flux_conv_inst;
     if ($sigma_inst>0.1) {
	 my $box=int(3*$sigma_inst*2/$dpix);
	 if ($box<3) {
	     $box=3;
	 }
	 my $kernel=zeroes(2*$box+1);
	 my    $norm=0;
	 for ($j=0;$j<2*$box+1;$j++) {
	     my $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	     set($kernel,$j,$gaus);
	     $norm=$norm+$gaus;
	 }
	 $kernel=$kernel/$norm;
	 $pdl_flux_conv_inst = conv2d $pdl_flux,$kernel;
     } else {
	 $pdl_flux_conv_inst=$pdl_flux;
     }

     
#
#
###################################################



#     print "After interp\n";
#     ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
#     print "end shift_convolve_pdl\n";
#
#    return $pdl_flux_l_wave_km_h;
    return $pdl_flux_conv_inst;

}


sub shift_convolve_pdl_NEW {
    my $pdl_wave_out=$_[0];
    my $pdl_wave=$_[1];
    my $pdl_flux=$_[2];
    my $redshift=$_[3];
    my $sigma_inst_in=$_[4];
    my $sigma_km_h=$_[5];
    my ($wmin,$wmax)=minmax(list($pdl_wave_out));
    my $dpix=$pdl_wave->at(1)-$pdl_wave->at(0);
    my $rsigma=$sigma_inst/$dpix;
    my $e=exp(1);
######################################################
# Convolved with Sigma_inst
#
#    my $box=int(3*$rsigma);
#
    my ($n_flux)=$pdl_flux->dims();
    my $pdl_ones=ones($n_flux);
    my $pdl_scale1=zeroes($n_flux);
    my $pdl_scale2=zeroes($n_flux);
    my $pdl_scale3=zeroes($n_flux);
    my $i;
    for ($i=0;$i<$n_flux;$i++) {
	my $val=1-$i/(0.5*$n_flux);
	if ($val>0) {
	    set($pdl_scale1,$i,$val);
	}
	my $val=$i/(0.5*$n_flux)-1;
	if ($val>0) {
	    set($pdl_scale2,$i,$val);
	}
    }
    $pdl_scale3=1-($pdl_scale2+$pdl_scale3);


# SIGMA_INST_0
#
    my $wmid=0.5*($wmax-$wmin);
     my $sigma_0=($sigma_km_h/300000)*$wmin;
     my $sigma_inst=sqrt($sigma_inst_in**2+$sigma_0**2);
    my $rsigma=$sigma_inst/$dpix;
     my $pdl_flux_conv_inst;
     my $box=int(5*$sigma_inst/$dpix);
     if ($box<3) {
	 $box=3;
     }
     my $kernel=zeroes(2*$box+1);
     my    $norm=0;
     for ($j=0;$j<2*$box+1;$j++) {
	 my $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	 set($kernel,$j,$gaus);
	 $norm=$norm+$gaus;
     }
     $kernel=$kernel/$norm;
     my $pdl_flux_conv_inst1 = conv2d $pdl_flux,$kernel;
     

# SIGMA_INST_1
#
     my $sigma_1=($sigma_km_h/300000)*$wmax;
     my $sigma_inst=sqrt($sigma_inst_in**2+$sigma_1**2);
   my $rsigma=$sigma_inst/$dpix;
#     my $pdl_flux_conv_inst;
#     my $box=int(3*3/$dpix);
     my $box=int(5*$sigma_inst/$dpix);
     if ($box<3) {
	 $box=3;
     }
     my $kernel=zeroes(2*$box+1);
     my    $norm=0;
     for ($j=0;$j<2*$box+1;$j++) {
	 my $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	 set($kernel,$j,$gaus);
	 $norm=$norm+$gaus;
     }
     $kernel=$kernel/$norm;
     my $pdl_flux_conv_inst2 = conv2d $pdl_flux,$kernel;

# SIGMA_INST_2
#
     my $sigma_2=($sigma_km_h/300000)*($wmax+$wmin)*0.5;
     my $sigma_inst=sqrt($sigma_inst_in**2+$sigma_2**2);
   my $rsigma=$sigma_inst/$dpix;
 #    my $pdl_flux_conv_inst;
#     my $box=int(3*3/$dpix);
     my $box=int(5*$sigma_inst/$dpix);
     if ($box<3) {
	 $box=3;
     }
     my $kernel=zeroes(2*$box+1);
     my    $norm=0;
     for ($j=0;$j<2*$box+1;$j++) {
	 my $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	 set($kernel,$j,$gaus);
	 $norm=$norm+$gaus;
     }
     $kernel=$kernel/$norm;
     my $pdl_flux_conv_inst3 = conv2d $pdl_flux,$kernel;

#
# Average mixing
#
#
#
#    print "SIGMA= $sigma_0,$sigma_1,$sigma_2 $sigma_inst\n";
    $pdl_flux_conv_inst1=$pdl_flux_conv_inst1*$pdl_scale1;
    $pdl_flux_conv_inst2=$pdl_flux_conv_inst2*$pdl_scale2;
    $pdl_flux_conv_inst3=$pdl_flux_conv_inst2*$pdl_scale3;
    my $pdl_flux_conv_inst=($pdl_flux_conv_inst1+$pdl_flux_conv_inst2+$pdl_flux_conv_inst3);

    my $pdl_flux_wave=interpol($pdl_wave_out,$pdl_wave,$pdl_flux_conv_inst);

    return $pdl_flux_wave;
    
}


sub print_time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    print "# TIME $sec $min $hour $mday $mon $year $wday $yday $isdst\n";
    print LOG "# TIME $sec $min $hour $mday $mon $year $wday $yday $isdst\n";
    my $sec_now=$hour*3600+$min*60+$sec;
    return $sec_now;
}

sub get_seconds {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my $sec_now=$yday*3600*24+$hour*3600+$min*60+$sec;
    return $sec_now;
}

sub FIT {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $k,$j,$i,$iter;

    my $iter_max=5;
    my $last_chi=1e12;
    
    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }


    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }
    $dpix_c_val[0]=$dpix_c_val[1];
$dpix_c=$wave_c[1]-$wave_c[0];

$rsigma=$sigma/$dpix_c;

    my @age_mod;
    my @met_mod;
for ($iii=0;$iii<$nf;$iii++) {
    $header="NAME".$iii;
    $name[$iii]=$pdl_flux_c_ini->hdr->{$header};
    $name_min=$name[$iii];
    $name_min =~ s/spec_ssp_//;
    $name_min =~ s/.spec//;    
    $name_min =~ s/.dat//;
    ($AGE,$MET)=split("_",$name_min);
    if ($AGE =~ "Myr") {
	$age=$AGE;
	$age =~ s/Myr//;
	$age=$age/1000;
    } else {
	$age=$AGE;
	$age =~ s/Gyr//;
    }
    $met=$MET;
    $met =~ s/z/0\./;

    $age_mod[$iii]=$age;
    $met_mod[$iii]=$met;
    $header="NORM".$iii;

 $val_ml=$pdl_flux_c_ini->hdr->{$header};
    if ($val_ml!=0) {
	$ml[$iii]=1/$val_ml;
    } else {
	$ml[$iii]=1;
    }
    
    $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    $kernel=zeroes(2*$box+1);
    $norm=0;
    $flux_c[$i]=0;
    for ($j=0;$j<2*$box+1;$j++) {
	$gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
    }
    $kernel=$kernel/$norm;
    

    $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;



    $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");


     my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
    ($n_c_out)=$out_spec_pdl->dims;


    for ($i=0;$i<$n_unc;$i++) {
	$val=$out_spec_pdl->at($i);
	if ($val eq "nan") {
	    $val=0;
	}
	$model[$i][$iii]=$val*$masked[$i];
	$model_no_mask[$i][$iii]=$val;#*$masked[$i];

	if ($masked[$i]>0) {
	    $error[$i]=0.01*abs($e_flux_unc[$i]);
	} else {
	    $error[$i]=0;
	}
	for ($j=0;$j<$nline;$j++) {
	    if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
#		    $error[$i]=0.1*$e_flux_unc[$i];
#		    $error[$i]=30;#*$e_flux_unc[$i];
		}
	}
    }
}

    $pdl_model=zeroes($n_unc,$nf);
    $pdl_model_no_mask=zeroes($n_unc,$nf);
    $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    $wave_res=$wave_unc[$i]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    $e_val=$error[$i];
	    $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
	}
    }
#
# We fit
#
    $pdl_flux_masked=pdl(@flux_masked);
    for ($j_mc=0;$j_mc<$n_mc;$j_mc++) {
	$C_left=1;
	$pdl_random=random($nf);
	# New converging criteria
#	$pdl_random=$pdl_C_last*

#	for ($j=0;$j<$nf;$j++) {
#	    $val=$pdl_random->at($j);
#	    $val=$val*$C_left;
#	    set($pdl_random,$j,$val);
#	    $C_left=$C_left-$val;	    
#	}
	$SUM=$pdl_random->sum;
	$pdl_random=$pdl_random/$SUM;
	$SUM=$pdl_random->sum;
	$pdl_C_now=zeroes($nf);

	$h_nf=1;
	for ($jj=0;$jj<$nf;$jj=$jj+$h_nf) {
	    $y_model_now=zeroes($n_unc);
	    $y_model_no_mask_now=zeroes($n_unc);
	    for ($j=0;$j<$nf;$j++) {
		$pdl_model_j=$pdl_model->slice(":,($j)");
		$pdl_model_no_mask_j=$pdl_model_no_mask->slice(":,($j)");
		$J=$j+$jj;
		if ($J>=$nf) {
		    $J=$J-$nf;
		}
		$val=$pdl_random->at($J);
		set($pdl_C_now,$j,$val);
		$y_model_now=$y_model_now+$val*$pdl_model_j;		
		$y_model_no_mask_now=$y_model_no_mask_now+$val*$pdl_model_no_mask_j;		
	    }
	    @dim_now=$y_model_now->dims;

	    $j1=int(0.4*$n_unc);
	    $j2=int(0.6*$n_unc);
	    my @a_norm;
	    my @b_norm;
	    my $j_a;
	    for ($j=$j1;$j<$j2;$j++) {
		$a_norm[$j_a]=$flux_unc[$j];
		$b_norm[$j_a]=$y_model_now->at($j);
		$j_a++;
	    }
	    $med_norm=median(@a_norm)/median(@b_norm);
	    $y_model_now=$y_model_now*$med_norm;
	    $y_model_no_mask_now=$y_model_no_mask_now*$med_norm;

##############################
# CHISQ VALUE
	    $chi=0;
	    $chi2=0;
	    $NFREE=0;
	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec[$j]=$y_model_now->at($j,0);
		$chi_sec[$j]=0;
		if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		    $chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/abs(1.5*$out_spec[$j]);
		    if ($have_error==0) {
			$chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($out_spec[$j]);
		    } else {
			$chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    }
		    $NFREE++;
		}
	    }
	    
	    $chi_sq=$chi;
	    if ($NFREE>0) {
		$chi_sq=($chi_sq/($NFREE))**0.5;
	    }

	    if ($chi_sq<$chi_sq_min_now) {
		$chi_sq_min_now=$chi_sq;
		$t=$coeffs->slice(":,0");
		$t .= $pdl_C_now;
		$y_model_end=$y_model_now;
		$y_model_no_mask_end=$y_model_no_mask_now;
	    }
	}
    }



    for ($j=0;$j<$n_unc;$j++) {
	$out_spec[$j]=$y_model_no_mask_end->at($j,0);
	$model_spec[$j]=$out_spec[$j];
	$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
	$model_spec_min[$j]=$model_spec[$j];
    }
    $chi_sq=$chi_sq_min_now;








    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
	print "\n";
	for ($j=0;$j<$n_unc;$j++) {
            $wave_res=$wave_unc[$j]/(1+$redshift);
            $dust_rat=A_l(3.1,$wave_res);
            $norm_C=0;
            $norm_C_mass=0;
            for ($k=0;$k<$nf;$k++) {
                $dust=10**(-0.4*$Av[$k]*$dust_rat);
                $C=$coeffs->at($k,0);
                $norm_C=$norm_C+$C;
                $norm_C_mass=$norm_C_mass+$C*$ml[$k];
            }
	}
	$age_min=0;
	$met_min=0;
	$Av_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	$Av_min_mass=0;
	open(C,">coeffs.out");
	print C "ID   AGE     MET    COEFF   Norm.Coeff  M/L   AV\n";
	print "---------------------------------------------\n";
	print "ID AGE     MET    COEFF   Norm.Coeff   M/L    AV\n";
	print "---------------------------------------------\n";

	$norm_C=0;
	$norm_C_mass=0;
	for ($k=0;$k<$nf;$k++) {
	    $dust=10**(-0.4*$Av[$k]*$dust_rat);
	    $C=$coeffs->at($k,0);
	    $norm_C=$norm_C+$C;
	    $norm_C_mass=$norm_C_mass+$C*$ml[$k];
	}


	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    if ($norm_C>0) {
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		$Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
		$CN=$C/$norm_C;
	    }
	    $C_now=$C*$med_norm;
	    printf(C "%-2d %-7.4f %-7.4f %-7.4f %-7.4f %-10.2f %-4.2f\n",$k,$age_mod[$k],$met_mod[$k],$C_now,$CN,$ml[$k],$Av[$k]);
	    if ($C>1e-5) {
		printf("%-2d %-7.4f %-7.4f %-7.4f %-7.4f %-10.2f %-4.2f\n",$k,$age_mod[$k],$met_mod[$k],$C_now,$CN,$ml[$k],$Av[$k]);
	    }
	    if ($norm_C_mass>0) {
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	    }

	}
	print "---------------------------------------------\n";
	close(C);
	$age_min=10**($age_min);
	$age_min_mass=10**($age_min_mass);
	
    }
    $name=$unc_file.", ";
    $scale="1";
    
    if ($plot==1) {

	pgbegin(0,$dev_plot,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($min_wave,$max_wave,$y_min,$y_max,0,0);
	pgsch(0.9);           # Set character height
	my $mean_Av=mean(@Av);
	pglabel("Wavelength","Counts","X=$chi_sq ($chi_joint) T=$age_min ($age_min_mass) Z=$met_min ($met_min_mass) Av=$mean_Av z=$redshift sigma=$sigma");
	pgsch(0.5);
	pgsci(1);
	pgline($n_unc,\@wave_unc,\@flux_unc);    
	pgsci(2);
	pgline($n_unc,\@wave_unc,\@out_spec);    
	pgsci(8);
	pgline($n_unc,\@wave_unc,\@model_spec_min);    
	pgsci(5);
	pgline($n_unc,\@wave_unc,\@res_spec);    

	if ($nc>0) {
	    pgsci(6);
	    pgline($nc,\@wave_clean,\@res_clean);    
	}
	pgsci(1);
	pgsci(13);
	pgsls(2);
	$j_plot=0;



	for ($j=0;$j<$nf;$j++) {
	    $C=$coeffs->at($j,0);
	    if ($C>0) {
		pgsci(13);
		my $model_now=$pdl_model->slice(",($j)");

		$model_now=$C*$model_now;
		@a_model_now=list($model_now);
		pgline($n_unc,\@wave_unc,\@a_model_now);    

	    }
	    if ($norm_C != 0) {
		$CN=$C/$norm_C;
	    } else {
		$CN=$C;
	    }
	    if ($C>0) {
		$pos_x=$min_wave*1.02;
		$pos_y=$y_max-0.05*($j_plot+1)*($y_max-$y_min);
		pgsch(0.8);           # Set character height
		if ($C>0) {
		    pgsci(1);
		} else {
		pgsci(15);
		}
		$CN=apr($CN);
		pgptxt($pos_x,$pos_y,0,0,"$CN $age_mod[$j] $met_mod[$j]");
		$j_plot++;
#	    print "$pos_x $pos_y\n"; <stdin>;
	    }
		pgsch(1.2);           # Set character height

	}
	pgsls(1);


	pgsci(1);


	pgclose;
	pgend;

    }

    return $chi_sq;
}


sub fit_ssp_rnd_guess_C {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
    my $pdl_C_input=$_[17];
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;
#    my $pdl
#    my $n_unc=$#flux_unc+1;
    my $n_unc=$#flux_unc+1;



    my $coeffs=zeroes($nf,3);    
    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);#$n_mc);

    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];

    my $rsigma=$sigma/$dpix_c;

    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	
	my $box=int(3*$rsigma);
	if ($box<3) {
	    $box=3;
	}
	my $kernel=zeroes(2*$box+1);
	my $norm=0;
	$flux_c[$i]=0;
	for ($j=0;$j<2*$box+1;$j++) {
	    $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	    set($kernel,$j,$gaus);
	    $norm=$norm+$gaus;
	}
	$kernel=$kernel/$norm;
        $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
	my $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	my ($n_c_out)=$out_spec_pdl->dims;
	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }


#    print "Av = @Av\n";


    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_good=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    my $pdl_dust=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    my $val=$model[$i][$j]*$dust;
 #	    print "$val $model[$i][$j]\n";
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
	    set($pdl_dust,$i,$j,$dust);
	}
#	print "Av[$j] = $Av[$j]\n";
    }




#
# We fit
#
    my $pdl_flux_masked=pdl(@flux_masked);
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);
    my $fact_q=0;
    for ($j_mc=0;$j_mc<$n_mc;$j_mc++) {
#	print "$j_mc/$n_mc\n";
	my $C_left=1;
	
#	if ($j_mc==0) {
#	    $fact_q=5;
#	}
#	if ($j_mc==1) {
#	    $fact_q=2;
#	}

	my $pdl_random=random($nf);
	my $pdl_grandom=random($nf);


	my $pdl_C_now;#=zeroes($nf);

	my $h_nf=1;



	my $jj=0;

	my $y_model_now=zeroes($n_unc);
	my $y_model_no_mask_now=zeroes($n_unc);
	my $sum_J=0;

	
	$pdl_C_now=$pdl_random->copy;

	for ($j=0;$j<$nf;$j++) {
	    my $val_random=2*$pdl_random->at($j)-1;
	    my $val_grandom=2*$pdl_grandom->at($j)-1;
	    my $C_val_input=$pdl_C_input->at($j);
	    my $v=$pdl_random;
#	    my $C_val=($C_val_input)+0.05*(2*$val_random-1)*$C_val_input+0.001*(2*$val_grandom-1);
#	    my $C_val=($C_val_input)+$fact_q*(0.2*$val_random*$C_val_input+0.1*$val_grandom);

#
	    my $C_val=($C_val_input)+$fact_q*(0.1*$val_grandom);
#	    my $C_val=($C_val_input)+(0.05*$val_grandom)+0.05*$fact_q/(0.5+$C_val_input);
	    
	 #   print "$j/$nf $v $val_random $val_grandom\n";

#	    my $C_val=($C_val_input)+0.1*(2*$val_grandom-1)/$nf;
	    if ($C_val>1) {
		$C_val=1;
	    }
	    if ($C_val<0) {
		$C_val=0;
	    }
	    set($pdl_C_now,$j,$C_val);
	    $sum_J=$sum_J+$C_val;
	}

#
	$pdl_C_now=$pdl_C_now/$sum_J;

#	print "PASO *\n";
	
	for ($j=0;$j<$nf;$j++) {
	    my $pdl_model_j=$pdl_model->slice(":,($j)");
	    my $pdl_model_no_mask_j=$pdl_model_no_mask->slice(":,($j)");	    
	    my $val=$pdl_C_now->at($j);#/$sum_J);
	    $y_model_now=$y_model_now+$val*$pdl_model_j;		
	    $y_model_no_mask_now=$y_model_no_mask_now+$val*$pdl_model_no_mask_j;		
	}
	my @dim_now=$y_model_now->dims;
	
#	print "PASO\n";

	my $j1=int(0.47*$n_unc);
	my $j2=int(0.53*$n_unc);
	
	my @a_norm;
	my @b_norm;
	my $j_a;
	for ($j=$j1;$j<$j2;$j++) {
	    $a_norm[$j_a]=$flux_unc[$j];
	    $b_norm[$j_a]=$y_model_now->at($j);
	    if (($a_norm[$j_a]>0)&&($b_norm[$j_a]>0)) {
		$j_a++;
	    }
	}
	
	my $med_b=median(@b_norm);
	my $med_norm;
	if ($med_b!=0) {
	    $med_norm=median(@a_norm)/median(@b_norm);
	} else {
	    $med_norm=1;
	}
	$MED_NORM[$jj]=$med_norm;
	$y_model_now=$y_model_now*$med_norm;
	$y_model_no_mask_now=$y_model_no_mask_now*$med_norm;
	#$pdl_C_now=$pdl_C_now/$med_norm;
	



##############################
# CHISQ VALUE
	my $chi=0;
	my $chi2=0;
	my $NFREE=0;
	my @out_spec;
	my @chi_sec;
	for ($j=0;$j<$n_unc;$j++) {
	    $out_spec[$j]=$y_model_now->at($j,0);
	    $chi_sec[$j]=0;
	    if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		$chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		if ($have_error==0) {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($out_spec[$j]);
		} else {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    }
		$NFREE++;
	    }
	}
	
	my $chi_sq=$chi;
	if ($NFREE>0) {
	    $chi_sq=($chi_sq/($NFREE))**0.5;
	}
	my $out_ps_now="junk";
	my $title="X = ".$chi_sq." Q=".$fact_q;



	plot_results($plot,pdl(@wave_unc),pdl(pdl(@flux_unc),pdl(@out_spec),pdl(@flux_unc)-pdl(@out_spec),$y_model_end),$out_ps_now,$title);
	print "$j_mc/$n_mc $chi_sq $chi_sq_min_now $val_random $val_grandom\n";
	#print "$chi_sq $pdl_C_now\n";

	    if ($chi_sq<$chi_sq_min_now) {

		$pdl_C_input=$pdl_C_now;
#	print "$pdl_C_input\n";
		$fact_q=0.9*$fact_q;
		if ($fact_q<0.1) {
		    $fact_q=0.3;
		}
		
#		print "INI_CAT $ini_cat/$n_mc\n";

		$chi_sq_min_now=$chi_sq;
		my $t=$coeffs->slice(":,0");
		$t .= $pdl_C_now;		
		$y_model_end=$y_model_now;
		$y_model_no_mask_end=$y_model_no_mask_now;
		my $nf_1=$nf-1;
		my $t=$coeffs_cat->slice("0:$nf_1,($ini_cat)");
		$t .= $pdl_C_now->slice(":,(0)");
		set($coeffs_cat,$nf,$ini_cat,$chi_sq);
		my $t=$pdl_model_spec_cat->slice(":,($ini_cat)");
		$t .= $y_model_no_mask_end;


		print "$coeffs\n";

		$ini_cat++;
	    }
#	    print "jj=$jj END\n";


# FOR JJ
#	}


#    print "PASO, $j_mc/$n_mc***\n";
    }

 #   print "Termino Bucle\n";

#    print "$j_mc/$n_mc\n";


#
# We construct the average model
# and the average coefficients
#

    my @model_spec;
    my @res_spec;
    my @model_spec_min;
    my @out_spec_now;
    my $SUM_W=0;
    my @out_coeffs;
    my @out_coeffs_e;
    my $N_CASES=0;
    for ($J=0;$J<$ini_cat;$J++) {
	my $CHI=$coeffs_cat->at($nf,$J);
	if ($CHI<1.1*$chi_sq_min_now) {
	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec_now[$j]=$out_spec_now[$j]+($pdl_model_spec_cat->at($j,$J))/$CHI;		
	    }	    
#	    my $pdl_coeff_now=$coeffs_cat->slice(":,$J");
	    for ($j=0;$j<$nf;$j++) {
		$out_coeffs[$j][$N_CASES]=$coeffs_cat->at($j,$J);#/$CHI;
	    }	    
	    $N_CASES++;
	    $SUM_W=$SUM_W+1/$CHI;
	}
    }




    if ($SUM_W==0) {
	$SUM_W=1;
    }
    for ($j=0;$j<$n_unc;$j++) {	
	$model_spec[$j]=$out_spec_now[$j]/$SUM_W;
	$out_spec[$j]=$out_spec_now[$j]/$SUM_W;
	$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
	$model_spec_min[$j]=$model_spec[$j];
    }

    $min_coeffs=$coeffs;

    
    for ($j=0;$j<$nf;$j++) {
	my @tmp;
	for ($J=0;$J<$N_CASES;$J++) {
	    $tmp[$J]=$out_coeffs[$j][$J];#/$SUM_W;
	}
	my $val=mean(@tmp);
	my $sigma=0.5*sigma(@tmp);
	my $sum_C=$pdl_C_input->sum;
	my $old_val=$coeffs->at($j,0);
#	my $old_val=$pdl_C_input->at($j)/$sum_C;#=$pdl_C_now->copy;
#	$sigma=$sigma/$sum_C
	set($coeffs,$j,0,$val);
	set($coeffs,$j,1,$sigma);
	set($coeffs,$j,2,$old_val);
    }


    my $chi_sq=$chi_sq_min_now;
	
    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);





    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
#	for ($j=0;$j<$n_unc;$j++) {
#            my $wave_res=$wave_unc[$j]/(1+$redshift);
#            my $dust_rat=A_l(3.1,$wave_res);
#            my $norm_C=0;
#            my $norm_C_mass=0;
#            for ($k=0;$k<$nf;$k++) {
#                my $dust=10**(-0.4*$Av[$k]*$dust_rat);
#                my $C=$coeffs->at($k,0);
#                $norm_C=$norm_C+$C;
#                $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#		set($coeffs_N,$k,0,$norm_C);
#		set($coeffs_NM,$k,0,$norm_C);
#            }
#	}


#	print "$chi_sq<$MIN_CHISQ\n $coeffs\n $coeffs_N\n";

	$age_min=0;
	$met_min=0;
	$Av_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	$Av_min_mass=0;

	my $norm_C=0;
	my $norm_C_mass=0;
	for ($k=0;$k<$nf;$k++) {
	    my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	    my $C=$coeffs->at($k,0);
	    $norm_C=$norm_C+$C;
	    $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#	    set($coeffs_N,$k,0,$norm_C);
#	    set($coeffs_NM,$k,0,$norm_C);
	}


	
	for ($k=0;$k<$nf;$k++) {
	    my $C=$coeffs->at($k,0);
	    set($coeffs_N,$k,0,$C/$norm_C);
	    set($coeffs_NM,$k,0,$C/$norm_C_mass);

	    if ($norm_C>0) {
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		$Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
		$CN=$C/$norm_C;
	    }
	    my $C_now=$C*$med_norm;
	    if ($norm_C_mass>0) {
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	    }

	}
	$age_min=10**($age_min);
	$age_min_mass=10**($age_min_mass);
	
    }
    my $name=$unc_file.", ";
    my $scale="1";
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);
    my $pdl_model_spec_min=pdl(@model_spec_min);
    my $pdl_res=pdl(@res_spec);

#    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);

    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,pdl(@wave_unc),pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);
    }


#    my $end_cat=$ini_cat-1;
#    my $coeffs_cat_sec=$coeffs_cat->slice(":,0:$end_cat");
    print "----- FINAL ----\n";
    print "$coeffs\n";
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res;
}


sub fit_ssp_rnd_guess_lin {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
#    my $pdl_C_input=$_[17];
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;
#    my $pdl
#    my $n_unc=$#flux_unc+1;
    my $n_unc=$#flux_unc+1;



    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);#$n_mc);

    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];
#    print "$dpix_c=$wave_c[1]-$wave_c[0]"; <stdin>;
    my $rsigma=$sigma/$dpix_c;

    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	
	my $box=int(3*$rsigma);
	if ($box<3) {
	    $box=3;
	}
	my $kernel=zeroes(2*$box+1);
	my $norm=0;
	$flux_c[$i]=0;
	for ($j=0;$j<2*$box+1;$j++) {
	    $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	    set($kernel,$j,$gaus);
	    $norm=$norm+$gaus;
	}
	$kernel=$kernel/$norm;
        $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
	my $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	my ($n_c_out)=$out_spec_pdl->dims;
	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }


#    print "Av = @Av\n";

    my $pdl_C_input=zeroes($nf);
    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_good=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    my $pdl_dust=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    my $val=$model[$i][$j]*$dust;
 #	    print "$val $model[$i][$j]\n";
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
#	    set($pdl_error,$i,(abs($val_now)**2));
	    set($pdl_dust,$i,$j,$dust);
	}
#	print "Av[$j] = $Av[$j]\n";
    }

    my $pdl_flux_masked=pdl(@flux_masked);

#######################################################
# LINEAR GUESS

# Just a linear FIT, without restrictions!
#    ($y_model_now, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;
    

    ($y_model_now, $coeffs) = linfit1d($pdl_flux_masked,$pdl_model,{Weights=>1/$pdl_error});


#
# We remove the models that are negative
#
    my $nf_new=0;
    $nf_neg=0;
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($C>0) {
	    $nf_new++;
	} else {
	    $nf_neg++;
	}
    }

    my @MOD;
    if ($nf_new>0) {
	while ($nf_neg>0) {
	    my $pdl_model_new=zeroes($n_unc,$nf_new);
	    my $nf_i=0;
	    for ($k=0;$k<$nf;$k++) {
		my $C=$coeffs->at($k,0);
		if ($C>0) {
		    my $t=$pdl_model_new->slice(":,$nf_i");
		    $t .= $pdl_model->slice(":,$k");
		    $MOD[$nf_i]=$k;
		    $nf_i++;
		} else {
		    set($coeffs,$k,0,0);
		}
	    }
	    my $coeffs_new;
	    my $yfit;
	    ($yfit, $coeffs_new) = linfit1d($pdl_flux_masked,$pdl_model_new,{Weights=>1/$pdl_error});

#    my $i_lin;    
#    my $n_lin=5;
#    for ($i_lin=0;$i_lin<$n_lin;$i_lin++) {
#	my $pdl_rnd=ones($n_unc)-2*random($n_unc);
#	my $pdl_flux_lin;
#	my $my_coeffs;
#	if ($i_lin==0) {
#	    $pdl_flux_lin=$pdl_flux_masked;
#	} else {
#	    $pdl_flux_lin=$pdl_flux_masked+$pdl_error*$pdl_rnd;
#	}
#	($y_model_now, $my_coeffs) = linfit1d($pdl_flux_lin,$pdl_model,{Weights=>1/$pdl_error});
#	print "$my_coeffs\n";
#	if ($i_lin==0) {
#	    $coeffs_new=$my_coeffs;
#	} else {
#	    $coeffs_new=$coeffs+$my_coeffs;
#	}
#
 #   }
 #   $coeffs_new=$coeffs_new/$n_lin;
##############################################






	    $y_model_now=$yfit;
	    my $nf_i=0;
	    $nf_neg=0;
	    $nf_new=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
		if ($C>0) {
		    my $val=$coeffs_new->at($nf_i,0);
		    $nf_i++;
		    if ($val>0) {
			set($coeffs,$k,0,$val);
			$nf_new++;
		    } else {
			set($coeffs,$k,0,0);
			$nf_neg++;
		    }
		}
	    }
		if ($nf_new==0) {
		    $nf_neg=0;
		}
	}
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    set($pdl_C_input,$k,$C);
	}

    } else {
	$nf_new=$nf;
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    set($pdl_C_input,$k,$C);
	}
    }

#
#
#



# End LINEAR GUESS
#######################################################
    
  #   $pdl_C_input=0.5*$pdl_C_input+0.5*random($nf);
     # $pdl_C_input=random($nf);
    my $sum_JUNK=$pdl_C_input->sum;
    $pdl_C_input=$pdl_C_input/$sum_JUNK;
    my $pdl_C_input_zero=$pdl_C_input->copy;

    my $min_C=1e12;
    my $max_C=0;

    for ($j=0;$j<$nf;$j++) {
	my $val_C=$pdl_C_input_zero->at($j);
	if ($val_C>0) {
	    if ($min_C<$val_C) {
		$min_C=$val_C;
	    }
	    if ($max_C>$val_C) {
		$max_C=$val_C;
	    }
	}
    }

#    print "$sum_JUNK $pdl_C_now \n";

    my $coeffs=zeroes($nf,3);    



#
# We fit
#
    my $ini_cat=0;
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);
    my $fact_q=1;
    my $i_mc=0;
    for ($j_mc=0;$j_mc<$n_mc;$j_mc++) {
#	print "$j_mc/$n_mc\n";
	my $C_left=1;
	
	if ($i_mc==0) {
	    $fact_q=0;
	} 
	if ($i_mc==1) {
	    $fact_q=1;
	} 

	$i_mc++;
#	if ($j_mc==1) {
#	    $fact_q=10;
#	}
#	srand(localtime);
	my $pdl_random=random($nf);
	my $pdl_grandom=random($nf);
	my $pdl_random_J=$nf*random($nf);

	my $pdl_C_now;#=zeroes($nf);

	my $h_nf=1;



	my $jj=0;

	my $y_model_now=zeroes($n_unc);
	my $y_model_no_mask_now=zeroes($n_unc);
	my $sum_J=0;

	
	$pdl_C_now=$pdl_random->copy;

	    if ($i_mc>1) {
		for ($j=0;$j<$nf;$j++) {
		    my $val_random=$pdl_random->at($j);
		    my $val_grandom=2*$pdl_grandom->at($j)-1;
		    my $C_val_input=$pdl_C_input->at($j);
		    my $C_val_zero=$pdl_C_input_zero->at($j);
		    
		    my $C_val=($C_val_input);
		    $C_val=$C_val+0.1*$fact_q*($val_grandom)/$nf;
		    if ($C_val>1) {
			$C_val=1;
		    }
		    if ($C_val<0) {
			$C_val=0;
		    }
		    set($pdl_C_now,$j,$C_val);
		}
	    } else {
#		$pdl_C_now=conv1d($pdl_C_input,pdl(0.01,0.05,0.9,0.05,0.01));
#		print "$pdl_C_input\n";
		$pdl_C_now=$pdl_C_input->copy;
		for ($j=0;$j<$nf;$j++) {
		    my $C_val_j=$pdl_C_input->at($j);
		    if ($C_val_j>0) {
			for ($jj=0;$jj<$nf;$jj++) {
			    if (($jj!=$j)&&($met_mod[$jj]==$met_mod[$j])) {
				my $C_val_jj=$pdl_C_input->at($jj);
#				print "(($jj!=$j)&&($met_mod[$jj]==$met_mod[$j])) $C_val_jj $C_val_j\n";
				$C_val_jj=$C_val_jj+0.1*$C_val_j/$nf;
				set($pdl_C_now,$jj,$C_val_jj);
				$C_val_j=$C_val_j-0.1*$C_val_j/$nf;
				set($pdl_C_now,$j,$C_val_j);
			    }
			}
		    }
		}
#		print "$pdl_C_now\n";
#		<stdin>;


	    }

#	print 
#
#

# We normalize!
#	

	my $sum_J=$pdl_C_now->sum;
	$pdl_C_now=$pdl_C_now/$sum_J;

#		my $sum_JUNK=$pdl_C_now->sum;
	#	print "$sum_JUNK $pdl_C_now \n";

#	print "PASO *\n";
	
	for ($j=0;$j<$nf;$j++) {
	    my $pdl_model_j=$pdl_model->slice(":,($j)");
	    my $pdl_model_no_mask_j=$pdl_model_no_mask->slice(":,($j)");	    
	    my $val=$pdl_C_now->at($j);#/$sum_J);
	    $y_model_now=$y_model_now+$val*$pdl_model_j;		
	    $y_model_no_mask_now=$y_model_no_mask_now+$val*$pdl_model_no_mask_j;		
	}
	my @dim_now=$y_model_now->dims;
	
#	print "PASO\n";

	my $j1=int(0.47*$n_unc);
	my $j2=int(0.53*$n_unc);
	
	my @a_norm;
	my @b_norm;
	my $j_a;
	for ($j=$j1;$j<$j2;$j++) {
	    $a_norm[$j_a]=$flux_unc[$j];
	    $b_norm[$j_a]=$y_model_now->at($j);
	    if (($a_norm[$j_a]>0)&&($b_norm[$j_a]>0)) {
		$j_a++;
	    }
	}
	
	my $med_b=median(@b_norm);
	my $med_norm;
	if ($med_b!=0) {
	    $med_norm=median(@a_norm)/median(@b_norm);
	} else {
	    $med_norm=1;
	}
	$MED_NORM[$jj]=$med_norm;
	$y_model_now=$y_model_now*$med_norm;
	$y_model_no_mask_now=$y_model_no_mask_now*$med_norm;
	#$pdl_C_now=$pdl_C_now/$med_norm;
	



##############################
# CHISQ VALUE
	my $chi=0;
	my $chi2=0;
	my $NFREE=0;
	my @out_spec;
	my @chi_sec;

	my $pdl_rand_noise=2*random($n_unc)-1;

	for ($j=0;$j<$n_unc;$j++) {
	    $out_spec[$j]=$y_model_now->at($j,0);
	    $chi_sec[$j]=0;
	    if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		my $rnd=0;#$e_flux_unc[$j]*$pdl_rand_noise->at($j);
		$chi=$chi+$masked[$j]*(($flux_masked[$j]+$rnd-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		if ($have_error==0) {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]+$rnd-$out_spec[$j])**2)/abs($out_spec[$j]);
		} else {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]+$rnd-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    }
		$NFREE++;
	    }
	}
	
	my $chi_sq=$chi;
	if ($NFREE>0) {
	    $chi_sq=($chi_sq/($NFREE))**0.5;
	}
	my $out_ps_now="junk";
	my $title="X = ".$chi_sq." Q=".$fact_q;



#	plot_results($plot,pdl(@wave_unc),pdl(pdl(@flux_unc),pdl(@out_spec),pdl(@flux_unc)-pdl(@out_spec),$y_model_end,pdl(@e_flux_unc)),$out_ps_now,$title);
#	print "$j_mc/$n_mc $chi_sq $chi_sq_min_now $fact_q\n";
	#print "$chi_sq $pdl_C_now\n";

	    if ($chi_sq<1.1*$chi_sq_min_now) {
#	    if ($chi_sq<1.15*$chi_sq_min_now) {
#		print "**** GOOD $chi_sq<1.1*$chi_sq_min_now)\n";
		$pdl_C_input=$pdl_C_now;
#	print "$pdl_C_input\n";
		$fact_q=0.95*$fact_q;
		$j_mc=0;
		if (($fact_q<0.05)&&($i_mc>1)) {
		    $j_mc=$n_mc;
		}
		
		if ($chi_sq<$chi_sq_min_now) {
		    $chi_sq_min_now=$chi_sq;
		}
		my $t=$coeffs->slice(":,0");
		$t .= $pdl_C_now;	
		my $sum_JUNK=$pdl_C_now->sum;
		$y_model_end=$y_model_now;
		$y_model_no_mask_end=$y_model_no_mask_now;
		my $nf_1=$nf-1;
		my $t=$coeffs_cat->slice("0:$nf_1,($ini_cat)");
		$t .= $pdl_C_now->slice(":,(0)");
		set($coeffs_cat,$nf,$ini_cat,$chi_sq);
		my $t=$pdl_model_spec_cat->slice(":,($ini_cat)");
		$t .= $y_model_no_mask_end;


		$ini_cat++;
		if ($ini_cat>$n_mc-2) {
		    $j_mc=$n_mc;
		}

	    }
#	    print "jj=$jj END\n";


# FOR JJ
#	}


#    print "PASO, $j_mc/$n_mc***\n";
    }

 #   print "Termino Bucle\n";

#    print "$j_mc/$n_mc\n";

#    print "End of Loop: INIT_CAT=$ini_cat\n";

#
# We construct the average model
# and the average coefficients
#

    my @model_spec;
    my @res_spec;
    my @model_spec_min;
    my @out_spec_now;
    my $SUM_W=0;
    my @out_coeffs;
    my @out_coeffs_e;
    my $N_CASES=0;
    my $pdl_model_final=zeroes($n_unc);
    for ($J=0;$J<$ini_cat;$J++) {
	my $CHI=$coeffs_cat->at($nf,$J);
	if ($CHI<1.1*$chi_sq_min_now) {

	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec_now[$j]=$out_spec_now[$j]+($pdl_model_spec_cat->at($j,$J))/$CHI;		
	    }

	    for ($j=0;$j<$nf;$j++) {
		my $val=$coeffs_cat->at($j,$J);#/$sum_J);
		$out_coeffs[$j][$N_CASES]=$val;#/$CHI;
#		my $pdl_model_no_mask_j=$pdl_model_no_mask->slice(":,($j)");	    
#		$pdl_model_final=$pdl_model_final+($val*$pdl_model_no_mask_j)/$CHI;
	    }	    
	    
	    $N_CASES++;
	    $SUM_W=$SUM_W+1/$CHI;
	}
    }

#    print "$N_CASES\n";

#		$out_spec_now[$j]=$out_spec_now[$j]+($pdl_model_spec_cat->at($j,$J))/$CHI;		



    if ($SUM_W==0) {
	$SUM_W=1;
    }
    for ($j=0;$j<$n_unc;$j++) {	
#	$out_spec_now[$j]=$pdl_model_final->at($j);
	$model_spec[$j]=$out_spec_now[$j]/$SUM_W;
	$out_spec[$j]=$out_spec_now[$j]/$SUM_W;
	$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
	$model_spec_min[$j]=$model_spec[$j];
    }

    $min_coeffs=$coeffs;

    
    for ($j=0;$j<$nf;$j++) {
	my @tmp;
	for ($J=0;$J<$N_CASES;$J++) {
	    $tmp[$J]=$out_coeffs[$j][$J];#/$SUM_W;
	}

	my $val=mean(@tmp);
	my $sigma=sigma(@tmp);
	my $sum_C=$pdl_C_input->sum;
#	print "J/NF = $j/$nf\n";
	my $old_val=$coeffs->at($j,0);
	set($coeffs,$j,0,$val);
	set($coeffs,$j,1,$sigma);
	set($coeffs,$j,2,$old_val);

    }



 #   print "PASO\n";

    my $chi_sq=$chi_sq_min_now;
	
#    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);





    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
#	for ($j=0;$j<$n_unc;$j++) {
#            my $wave_res=$wave_unc[$j]/(1+$redshift);
#            my $dust_rat=A_l(3.1,$wave_res);
#            my $norm_C=0;
#            my $norm_C_mass=0;
#            for ($k=0;$k<$nf;$k++) {
#                my $dust=10**(-0.4*$Av[$k]*$dust_rat);
#                my $C=$coeffs->at($k,0);
#                $norm_C=$norm_C+$C;
#                $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#		set($coeffs_N,$k,0,$norm_C);
#		set($coeffs_NM,$k,0,$norm_C);
#            }
#	}


#	print "$chi_sq<$MIN_CHISQ\n $coeffs\n $coeffs_N\n";

	$age_min=0;
	$met_min=0;
	$Av_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	$Av_min_mass=0;

	my $norm_C=0;
	my $norm_C_mass=0;
	for ($k=0;$k<$nf;$k++) {
	    my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	    my $C=$coeffs->at($k,0);
	    $norm_C=$norm_C+$C;
	    $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#	    set($coeffs_N,$k,0,$norm_C);
#	    set($coeffs_NM,$k,0,$norm_C);
	}


	
	for ($k=0;$k<$nf;$k++) {
	    my $C=$coeffs->at($k,0);
	    set($coeffs_N,$k,0,$C/$norm_C);
	    set($coeffs_NM,$k,0,$C/$norm_C_mass);

	    if ($norm_C>0) {
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		$Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
		$CN=$C/$norm_C;
	    }
	    my $C_now=$C*$med_norm;
	    if ($norm_C_mass>0) {
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	    }

	}
	$age_min=10**($age_min);
	$age_min_mass=10**($age_min_mass);
	
    }
    my $name=$unc_file.", ";
    my $scale="1";
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);
    my $pdl_model_spec_min=pdl(@model_spec_min);
    my $pdl_res=pdl(@res_spec);

#    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);

    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,pdl(@wave_unc),pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);
    }


#    my $end_cat=$ini_cat-1;
#    my $coeffs_cat_sec=$coeffs_cat->slice(":,0:$end_cat");
#    print "----- FINAL ----\n";
#    print "$coeffs\n";
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res,$pdl_C_input_zero;
}


sub fit_ssp_lin_MC {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
    my $sigma_inst=$_[17];
#    my $pdl_C_input=$_[17];
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;
#    my $pdl
#    my $n_unc=$#flux_unc+1;
    my $n_unc=$#flux_unc+1;



    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);#$n_mc);

    my $pdl_1st_model;

    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];
#    print "$dpix_c=$wave_c[1]-$wave_c[0]"; <stdin>;


    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	

	my $pdl_flux_c_ini_now=$pdl_flux_c_ini->slice(",($iii)");
	my $out_spec_pdl=shift_convolve_pdl(pdl(@wave_unc),pdl(@wave_c),$pdl_flux_c_ini_now,0,$sigma_inst,$sigma);
#        $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
#	my $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
#	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	my ($n_c_out)=$out_spec_pdl->dims;
	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }


#    print "Av = @Av\n";

    my $pdl_C_input=zeroes($nf);
    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_good=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    my $pdl_dust=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    my $val=$model[$i][$j]*$dust;
 #	    print "$val $model[$i][$j]\n";
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
#	    set($pdl_error,$i,(abs($val_now)**2));
	    set($pdl_dust,$i,$j,$dust);
	}
#	print "Av[$j] = $Av[$j]\n";
    }

    my $pdl_flux_masked=pdl(@flux_masked);

#######################################################
# LINEAR GUESS

    my $j_iter=0;
    my $n_iter=$n_mc;
    my $coeffs_iter=zeroes($nf,$n_iter);
    my $pdl_C_rms;

    for ($j_iter=0;$j_iter<$n_iter;$j_iter++) {
	my $pdl_gr=grandom($n_unc);
	$pdl_gr->inplace->clip(-1,1);
	my $pdl_noise=(1/$pdl_error)*$pdl_gr;
	my $pdl_flux_to_fit=$pdl_flux_masked+$pdl_noise;
	($y_model_now, $coeffs) = linfit1d($pdl_flux_to_fit,$pdl_model,{Weights=>$pdl_error});

	$pdl_1st_model=$y_model_now->copy();

#
# We remove the models that are negative
#
    my $nf_new=0;
    $nf_neg=0;
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($C>0) {
	    $nf_new++;
	} else {
	    $nf_neg++;
	}
    }

    my @MOD;
    if ($nf_new>0) {
	while ($nf_neg>0) {
	    my $pdl_model_new=zeroes($n_unc,$nf_new);
	    my $nf_i=0;
	    for ($k=0;$k<$nf;$k++) {
		my $C=$coeffs->at($k,0);
		if ($C>0) {
		    my $t=$pdl_model_new->slice(":,$nf_i");
		    $t .= $pdl_model->slice(":,$k");
		    $MOD[$nf_i]=$k;
		    $nf_i++;
		} else {
		    set($coeffs,$k,0,0);
		}
	    }
	    my $coeffs_new;
	    my $yfit;
	    ($yfit, $coeffs_new) = linfit1d($pdl_flux_to_fit,$pdl_model_new,{Weights=>$pdl_error});

	    $pdl_1st_model=$yfit->copy();

	    $y_model_now=$yfit;
	    my $nf_i=0;
	    $nf_neg=0;
	    $nf_new=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
		if ($C>0) {
		    my $val=$coeffs_new->at($nf_i,0);
		    $nf_i++;
		    if ($val>0) {
			set($coeffs,$k,0,$val);
			$nf_new++;
		    } else {
			set($coeffs,$k,0,0);
			$nf_neg++;
		    }
		}
	    }
		if ($nf_new==0) {
		    $nf_neg=0;
		}
	}
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
#	    set($pdl_C_input,$k,$C);
	}

    } else {
	$nf_new=$nf;
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
#	    set($pdl_C_input,$k,$C);
	}
    }

#
#
#
#	print "J_iter/N_iter = $j_iter/$n_iter\n";

	my $t=$coeffs_iter->slice(":,$j_iter");
	$t .= $coeffs->slice(":,0");

#
#	my $out_ps_now="junk.ps";
#	my $title="$j_iter\n";
#	plot_results($plot,pdl(@wave_unc),pdl($pdl_flux_to_fit,$pdl_flux_masked,$y_model_now,$pdl_flux_to_fit-$y_model_now),$out_ps_now,$title);
#	print "$t\n";
		     
#	print "Press Enter"; <stdin>;

		   
#    my $coeffs_iter=zeroes($nf,$n_iter);
    }

#    print "$coeffs_iter\n";

# End LINEAR GUESS
#######################################################
    
  #   $pdl_C_input=0.5*$pdl_C_input+0.5*random($nf);
     # $pdl_C_input=random($nf);

    $pdl_C_input = sumover $coeffs_iter->xchg(0,1);
    my $sum_JUNK=$pdl_C_input->sum;
    $pdl_C_input=$pdl_C_input/$sum_JUNK;
    my $pdl_C_input_zero=$pdl_C_input->copy;
    $coeffs_iter=$coeffs_iter/$sum_JUNK;


    $pdl_C_rms=$pdl_C_input->copy;
    for ($j=0;$j<$nf;$j++) {
	my $sec=$coeffs_iter->slice("($j),:");
	my @stats_sec=stats($sec);
	set($pdl_C_rms,$j,$stats_sec[1]);
    }

#    print "$pdl_C_input\n $pdl_C_rms\n"; <stdin>;


    my $min_C=1e12;
    my $max_C=0;

    for ($j=0;$j<$nf;$j++) {
	my $val_C=$pdl_C_input_zero->at($j);
	if ($val_C>0) {
	    if ($min_C<$val_C) {
		$min_C=$val_C;
	    }
	    if ($max_C>$val_C) {
		$max_C=$val_C;
	    }
	}
    }

#    print "$sum_JUNK $pdl_C_now \n";

    my $coeffs=zeroes($nf,3);    



#
# We fit
#
    my $ini_cat=0;
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);
    my $fact_q=1;
    my $i_mc=0;
    for ($j_mc=0;$j_mc<$n_mc;$j_mc++) {
#	print "$j_mc/$n_mc\n";
	my $C_left=1;
	
	if ($i_mc==0) {
	    $fact_q=0;
	} 
	if ($i_mc==1) {
	    $fact_q=1;
	} 

	$i_mc++;
#	if ($j_mc==1) {
#	    $fact_q=10;
#	}
#	srand(localtime);
	my $pdl_random=random($nf);
	my $pdl_grandom=random($nf);
	my $pdl_random_J=$nf*random($nf);

	my $pdl_C_now;#=zeroes($nf);

	my $h_nf=1;



	my $jj=0;

	my $y_model_now=zeroes($n_unc);
	my $y_model_no_mask_now=zeroes($n_unc);
	my $sum_J=0;

	
	$pdl_C_now=$pdl_random->copy;
	
	if ($i_mc>1) {
	    for ($j=0;$j<$nf;$j++) {
		my $val_random=$pdl_random->at($j);
		my $val_grandom=2*$pdl_grandom->at($j)-1;
		my $C_val_input=$pdl_C_input->at($j);
		my $C_val_zero=$pdl_C_input_zero->at($j);
		
		my $C_val=$C_val_zero;
#		my $C_val=$C_val_input;
		my $C_rms=$pdl_C_rms->at($j);
		$C_val=$C_val+2*$fact_q*$C_rms;
		if ($C_val>1) {
		    $C_val=1;
		}
		if ($C_val<0) {
		    $C_val=0;
		}
		set($pdl_C_now,$j,$C_val);
	    }
	} else {
	    $pdl_C_now=$pdl_C_input->copy;
	}
	
#	print 
#
#

# We normalize!
#	

	my $sum_J=$pdl_C_now->sum;
	$pdl_C_now=$pdl_C_now/$sum_J;

#		my $sum_JUNK=$pdl_C_now->sum;
	#	print "$sum_JUNK $pdl_C_now \n";

#	print "PASO *\n";
	
	for ($j=0;$j<$nf;$j++) {
	    my $pdl_model_j=$pdl_model->slice(":,($j)");
	    my $pdl_model_no_mask_j=$pdl_model_no_mask->slice(":,($j)");	    
	    my $val=$pdl_C_now->at($j);#/$sum_J);
	    $y_model_now=$y_model_now+$val*$pdl_model_j;		
	    $y_model_no_mask_now=$y_model_no_mask_now+$val*$pdl_model_no_mask_j;		
	}
	my @dim_now=$y_model_now->dims;
	
#	print "PASO\n";

	my $j1=int(0.47*$n_unc);
	my $j2=int(0.53*$n_unc);
	
	my @a_norm;
	my @b_norm;
	my $j_a;
	for ($j=$j1;$j<$j2;$j++) {
	    $a_norm[$j_a]=$flux_unc[$j];
	    $b_norm[$j_a]=$y_model_now->at($j);
	    if (($a_norm[$j_a]>0)&&($b_norm[$j_a]>0)) {
		$j_a++;
	    }
	}
	
	my $med_b=median(@b_norm);
	my $med_norm;
	if ($med_b!=0) {
	    $med_norm=median(@a_norm)/median(@b_norm);
	} else {
	    $med_norm=1;
	}
	$MED_NORM[$jj]=$med_norm;
	$y_model_now=$y_model_now*$med_norm;
	$y_model_no_mask_now=$y_model_no_mask_now*$med_norm;
	#$pdl_C_now=$pdl_C_now/$med_norm;
	



##############################
# CHISQ VALUE
	my $chi=0;
	my $chi2=0;
	my $NFREE=0;
	my @out_spec;
	my @chi_sec;

	my $pdl_rand_noise=2*random($n_unc)-1;

	
#	my $smooth_ratio=smooth_ratio(pdl(@flux_unc),$y_model_now->slice(":,0"),int($sigma));
	

	for ($j=0;$j<$n_unc;$j++) {
#	    $out_spec[$j]=($y_model_now->at($j,0))*($smooth_ratio->at($j));
	    $out_spec[$j]=$y_model_now->at($j,0);
	    $chi_sec[$j]=0;
	    if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		my $rnd=0;#$e_flux_unc[$j]*$pdl_rand_noise->at($j);
		$chi=$chi+$masked[$j]*(($flux_masked[$j]+$rnd-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		if ($have_error==0) {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]+$rnd-$out_spec[$j])**2)/abs($out_spec[$j]);
		} else {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]+$rnd-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    }
		$NFREE++;
	    }
	}
	
	my $chi_sq=$chi;
	if ($NFREE>0) {
	    $chi_sq=($chi_sq/($NFREE))**0.5;
	}
	my $out_ps_now="junk";
	my $title="X = ".$chi_sq." Q=".$fact_q;


	print "CHI_SQ_LOOP = $chi_sq\n";

	plot_results($plot,pdl(@wave_unc),pdl(pdl(@flux_unc),pdl(@out_spec),pdl(@flux_unc)-pdl(@out_spec),$y_model_end,pdl(@e_flux_unc)),$out_ps_now,$title);

	    if ($chi_sq<1.1*$chi_sq_min_now) {
		$pdl_C_input=$pdl_C_now;
		$fact_q=0.95*$fact_q;
		$j_mc=0;
		if (($fact_q<0.05)&&($i_mc>1)) {
		    $j_mc=$n_mc;
		}
		
		if ($chi_sq<$chi_sq_min_now) {
		    $chi_sq_min_now=$chi_sq;
		}
		my $t=$coeffs->slice(":,0");
		$t .= $pdl_C_now;	
		my $sum_JUNK=$pdl_C_now->sum;
		$y_model_end=$y_model_now;
		$y_model_no_mask_end=$y_model_no_mask_now;
		my $nf_1=$nf-1;
		my $t=$coeffs_cat->slice("0:$nf_1,($ini_cat)");
		$t .= $pdl_C_now->slice(":,(0)");
		set($coeffs_cat,$nf,$ini_cat,$chi_sq);
		my $t=$pdl_model_spec_cat->slice(":,($ini_cat)");
		$t .= $y_model_no_mask_end;


		$ini_cat++;
		if ($ini_cat>$n_mc-2) {
		    $j_mc=$n_mc;
		}

	    }
#	    print "jj=$jj END\n";


# FOR JJ
#	}


#    print "PASO, $j_mc/$n_mc***\n";
    }

 #   print "Termino Bucle\n";

#    print "$j_mc/$n_mc\n";

#    print "End of Loop: INIT_CAT=$ini_cat\n";

#
# We construct the average model
# and the average coefficients
#

    my @model_spec;
    my @res_spec;
    my @model_spec_min;
    my @out_spec_now;
    my $SUM_W=0;
    my @out_coeffs;
    my @out_coeffs_e;
    my $N_CASES=0;
    my $pdl_model_final=zeroes($n_unc);
    for ($J=0;$J<$ini_cat;$J++) {
	my $CHI=$coeffs_cat->at($nf,$J);
	if ($CHI<1.1*$chi_sq_min_now) {

	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec_now[$j]=$out_spec_now[$j]+($pdl_model_spec_cat->at($j,$J))/$CHI;		
	    }

	    for ($j=0;$j<$nf;$j++) {
		my $val=$coeffs_cat->at($j,$J);#/$sum_J);
		$out_coeffs[$j][$N_CASES]=$val;#/$CHI;
	    }	    
	    
	    $N_CASES++;
	    $SUM_W=$SUM_W+1/$CHI;
	}
    }

#    print "$N_CASES\n";

#		$out_spec_now[$j]=$out_spec_now[$j]+($pdl_model_spec_cat->at($j,$J))/$CHI;		



    if ($SUM_W==0) {
	$SUM_W=1;
#
# No better solution found than the 1st one!!!
#
	for ($j=0;$j<$n_unc;$j++) {
	    my $val=$pdl_1st_model->at($j);
	    $model_spec[$j]=$val;
	    $out_spec[$j]=$val;
	    $res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
	    $model_spec_min[$j]=$model_spec[$j];
	}
	
    } else {
	for ($j=0;$j<$n_unc;$j++) {	
	    $model_spec[$j]=$out_spec_now[$j]/$SUM_W;
	    $out_spec[$j]=$out_spec_now[$j]/$SUM_W;
	    $res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
	    $model_spec_min[$j]=$model_spec[$j];
	}
    }

    $min_coeffs=$coeffs;

    
    for ($j=0;$j<$nf;$j++) {
	my @tmp;
	for ($J=0;$J<$N_CASES;$J++) {
	    $tmp[$J]=$out_coeffs[$j][$J];#/$SUM_W;
	}

	my $val=mean(@tmp);
	my $sigma=sigma(@tmp);       
	my $sigma_MC=$pdl_C_rms->at($j);
	$sigma=sqrt($sigma**2+$sigma_MC**2);
#	print "$sigma $sigma_MC\n";
	my $sum_C=$pdl_C_input->sum;
#	print "J/NF = $j/$nf\n";
	my $old_val=$coeffs->at($j,0);
	set($coeffs,$j,0,$val);
	set($coeffs,$j,1,$sigma);
	set($coeffs,$j,2,$old_val);

    }



 #   print "PASO\n";

    my $chi_sq=$chi_sq_min_now;
	
#    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);





    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
#	for ($j=0;$j<$n_unc;$j++) {
#            my $wave_res=$wave_unc[$j]/(1+$redshift);
#            my $dust_rat=A_l(3.1,$wave_res);
#            my $norm_C=0;
#            my $norm_C_mass=0;
#            for ($k=0;$k<$nf;$k++) {
#                my $dust=10**(-0.4*$Av[$k]*$dust_rat);
#                my $C=$coeffs->at($k,0);
#                $norm_C=$norm_C+$C;
#                $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#		set($coeffs_N,$k,0,$norm_C);
#		set($coeffs_NM,$k,0,$norm_C);
#            }
#	}


#	print "$chi_sq<$MIN_CHISQ\n $coeffs\n $coeffs_N\n";

	$age_min=0;
	$met_min=0;
	$Av_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	$Av_min_mass=0;

	my $norm_C=0;
	my $norm_C_mass=0;
	for ($k=0;$k<$nf;$k++) {
	    my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	    my $C=$coeffs->at($k,0);
	    $norm_C=$norm_C+$C;
	    $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#	    set($coeffs_N,$k,0,$norm_C);
#	    set($coeffs_NM,$k,0,$norm_C);
	}


	
	for ($k=0;$k<$nf;$k++) {
	    my $C=$coeffs->at($k,0);
	    set($coeffs_N,$k,0,$C/$norm_C);
	    set($coeffs_NM,$k,0,$C/$norm_C_mass);

	    if ($norm_C>0) {
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		$Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
		$CN=$C/$norm_C;
	    }
	    my $C_now=$C*$med_norm;
	    if ($norm_C_mass>0) {
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	    }

	}
	$age_min=10**($age_min);
	$age_min_mass=10**($age_min_mass);
	
    }
    my $name=$unc_file.", ";
    my $scale="1";
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);
    my $pdl_model_spec_min=pdl(@model_spec_min);
#    my $smooth_ratio=smooth_ratio(pdl(@flux_unc),$y_model_spec_min,int($sigma));
#    $pdl_model_spec_min=$pdl_model_spec_min*$smooth_ratio;

    my $pdl_res=pdl(@res_spec);

#    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);

    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,pdl(@wave_unc),pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);
    }


#    my $end_cat=$ini_cat-1;
#    my $coeffs_cat_sec=$coeffs_cat->slice(":,0:$end_cat");
#    print "----- FINAL ----\n";
#    print "$coeffs\n";
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res,$pdl_C_input_zero,$pdl_C_rms; 
}


sub fit_ssp_rnd {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
    my $pdl_rat_now=$_[17];
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;

#    my $n_unc=$#flux_unc+1;
    my $n_unc=$#flux_unc+1;



    my $coeffs=zeroes($nf,3);    
    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc*$n_mc);

    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];

    my $rsigma=$sigma/$dpix_c;

    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	
	my $box=int(3*$rsigma);
	if ($box<3) {
	    $box=3;
	}
	my $kernel=zeroes(2*$box+1);
	my $norm=0;
	$flux_c[$i]=0;
	for ($j=0;$j<2*$box+1;$j++) {
	    $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	    set($kernel,$j,$gaus);
	    $norm=$norm+$gaus;
	}
	$kernel=$kernel/$norm;
        $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
	my $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	my ($n_c_out)=$out_spec_pdl->dims;
	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }

#    print "Av = @Av\n";


    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    my $pdl_dust=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    my $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
	    set($pdl_dust,$i,$j,$dust);
	}
#	print "Av[$j] = $Av[$j]\n";
    }
#
# We fit
#


# 
   my $pdl_flux_masked=pdl(@flux_masked);
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);
    for ($j_mc=0;$j_mc<$n_mc;$j_mc++) {
	my $C_left=1;

	my $pdl_random=random($nf);
	for ($j=0;$j<$nf;$j++) {
	    my $val=$pdl_random->at($j);
#	    if ($val<0.01) {
#		$val=0;
#	    }
	    my $val=$val*$C_left;
	    set($pdl_random,$j,$val);
	    $C_left=$C_left-$val;	    
	}
	my $SUM=$pdl_random->sum;
	my $pdl_random=$pdl_random/$SUM;
	my $SUM=$pdl_random->sum;
	my $pdl_C_now;#=zeroes($nf);

	my $h_nf=1;




	for ($jj=0;$jj<$nf;$jj=$jj+$h_nf) {
	    my $y_model_now=zeroes($n_unc);
	    my $y_model_no_mask_now=zeroes($n_unc);
	    my $sum_J=0;
	    $pdl_C_now=$pdl_random->copy;
	    for ($j=0;$j<$nf;$j++) {
	#	my $J=$j+$jj;
	#	if ($J>=$nf) {
	#	    $J=$J-$nf;
	#	}
		my $pdl_J=$nf*random($nf);
		my $J=floor($pdl_J->at($j));
		my $val_J=$pdl_C_now->at($J);
		my $val_j=$pdl_C_now->at($j);
		set($pdl_C_now,$j,$val_J);
		set($pdl_C_now,$J,$val_j);
#		$sum_J=$sum_J+$val;
		#
		# We center in a good solution
		#
	    }

	    for ($j=0;$j<$nf;$j++) {
		my $pdl_model_j=$pdl_model->slice(":,($j)");
		my $pdl_model_no_mask_j=$pdl_model_no_mask->slice(":,($j)");

		my $val=$pdl_C_now->at($j);#/$sum_J);
		$y_model_now=$y_model_now+$val*$pdl_model_j;		
		$y_model_no_mask_now=$y_model_no_mask_now+$val*$pdl_model_no_mask_j;		
	    }
	    my @dim_now=$y_model_now->dims;

	    my $j1=int(0.47*$n_unc);
	    my $j2=int(0.53*$n_unc);

#	    my $j1=int(0.17*$n_unc);
#	    my $j2=int(0.23*$n_unc);

#	    my $j1=int(0.67*$n_unc);
#	    my $j2=int(0.73*$n_unc);

#	    my $j1=0;#int(0.4*$n_unc);
#	    my $j2=$n_unc;#int(0.6*$n_unc);
	    my @a_norm;
	    my @b_norm;
	    my $j_a;
	    for ($j=$j1;$j<$j2;$j++) {
		$a_norm[$j_a]=$flux_unc[$j];
		$b_norm[$j_a]=$y_model_now->at($j);
		if (($a_norm[$j_a]>0)&&($b_norm[$j_a]>0)) {
		    $j_a++;
		}
	    }

	    my $med_b=median(@b_norm);
#	    my $med_b=mean(@b_norm);
	    my $med_norm;
	    if ($med_b!=0) {
		$med_norm=median(@a_norm)/median(@b_norm);
#		$med_norm=mean(@a_norm)/mean(@b_norm);
	    } else {
		$med_norm=1;
	    }
#	    my $med_norm=mean(@a_norm)/mean(@b_norm);
	    $MED_NORM[$jj]=$med_norm;
	    $y_model_now=$y_model_now*$med_norm;
	    $y_model_no_mask_now=$y_model_no_mask_now*$med_norm;

##############################
# CHISQ VALUE
	    my $chi=0;
	    my $chi2=0;
	    my $NFREE=0;
	    my @out_spec;
	    my @chi_sec;
	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec[$j]=$y_model_now->at($j,0);
		$chi_sec[$j]=0;
		if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		    $chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    if ($have_error==0) {
			$chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($out_spec[$j]);
		    } else {
			$chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    }
		    $NFREE++;
		}
	    }
	    
	    my $chi_sq=$chi;
	    if ($NFREE>0) {
		$chi_sq=($chi_sq/($NFREE))**0.5;
	    }

#	    print "jj=$jj\n";
#	    print "$jj $j_mc/$n_mc $chi_sq<$chi_sq_min_now $NFREE $n_unc\n";
	    if ($chi_sq<$chi_sq_min_now) {
#		print "INI_CAT $ini_cat/$n_mc\n";

		$chi_sq_min_now=$chi_sq;
		my $t=$coeffs->slice(":,0");
		$t .= $pdl_C_now;
		$y_model_end=$y_model_now;
		$y_model_no_mask_end=$y_model_no_mask_now;
		my $nf_1=$nf-1;
		my $t=$coeffs_cat->slice("0:$nf_1,($ini_cat)");
		$t .= $pdl_C_now->slice(":,(0)");
		set($coeffs_cat,$nf,$ini_cat,$chi_sq);
		my $t=$pdl_model_spec_cat->slice(":,($ini_cat)");
		$t .= $y_model_no_mask_end;

		$ini_cat++;
	    }
#	    print "jj=$jj END\n";
	}
#    print "PASO, $j_mc/$n_mc***\n";
    }

 #   print "Termino Bucle\n";

#    print "$j_mc/$n_mc\n";


#
# We construct the average model
# and the average coefficients
#

    my @model_spec;
    my @res_spec;
    my @model_spec_min;
    my @out_spec_now;
    my $SUM_W=0;
    my @out_coeffs;
    my @out_coeffs_e;
    my $N_CASES=0;
    for ($J=0;$J<$ini_cat;$J++) {
	my $CHI=$coeffs_cat->at($nf,$J);
	if ($CHI<2*$chi_sq_min_now) {
	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec_now[$j]=$out_spec_now[$j]+($pdl_model_spec_cat->at($j,$J))/$CHI;		
	    }
#	    my $pdl_coeff_now=$coeffs_cat->slice(":,$J");
	    for ($j=0;$j<$nf;$j++) {
		$out_coeffs[$j][$N_CASES]=$coeffs_cat->at($j,$J);#/$CHI;
	    }
	    $N_CASES++;
	    $SUM_W=$SUM_W+1/$CHI;
	}
    }

    for ($j=0;$j<$n_unc;$j++) {
	$model_spec[$j]=$out_spec_now[$j]/$SUM_W;
	$out_spec[$j]=$out_spec_now[$j]/$SUM_W;
	$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
	$model_spec_min[$j]=$model_spec[$j];
    }

    $min_coeffs=$coeffs;

    for ($j=0;$j<$nf;$j++) {
	my @tmp;
	for ($J=0;$J<$N_CASES;$J++) {
	    $tmp[$J]=$out_coeffs[$j][$J];#/$SUM_W;
#	    print "$j,$J,$tmp[$J]\n";
	}
	my $val=mean(@tmp);
	my $sigma=0.5*sigma(@tmp);
	my $old_val=$coeffs->at($j,0);
	set($coeffs,$j,0,$val);
	set($coeffs,$j,1,$sigma);
	set($coeffs,$j,2,$old_val);
    }

    my $chi_sq=$chi_sq_min_now;
	
    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);





    if ($chi_sq<$MIN_CHISQ) {
	$MIN_CHISQ=$chi_sq;
#	for ($j=0;$j<$n_unc;$j++) {
#            my $wave_res=$wave_unc[$j]/(1+$redshift);
#            my $dust_rat=A_l(3.1,$wave_res);
#            my $norm_C=0;
#            my $norm_C_mass=0;
#            for ($k=0;$k<$nf;$k++) {
#                my $dust=10**(-0.4*$Av[$k]*$dust_rat);
#                my $C=$coeffs->at($k,0);
#                $norm_C=$norm_C+$C;
#                $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#		set($coeffs_N,$k,0,$norm_C);
#		set($coeffs_NM,$k,0,$norm_C);
#            }
#	}


#	print "$chi_sq<$MIN_CHISQ\n $coeffs\n $coeffs_N\n";

	$age_min=0;
	$met_min=0;
	$Av_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	$Av_min_mass=0;

	my $norm_C=0;
	my $norm_C_mass=0;
	for ($k=0;$k<$nf;$k++) {
	    my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	    my $C=$coeffs->at($k,0);
	    $norm_C=$norm_C+$C;
	    $norm_C_mass=$norm_C_mass+$C*$ml[$k];
#	    set($coeffs_N,$k,0,$norm_C);
#	    set($coeffs_NM,$k,0,$norm_C);
	}


	
	for ($k=0;$k<$nf;$k++) {
	    my $C=$coeffs->at($k,0);
	    set($coeffs_N,$k,0,$C/$norm_C);
	    set($coeffs_NM,$k,0,$C/$norm_C_mass);

	    if ($norm_C>0) {
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		$Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
		$CN=$C/$norm_C;
	    }
	    my $C_now=$C*$med_norm;
	    if ($norm_C_mass>0) {
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	    }

	}
	$age_min=10**($age_min);
	$age_min_mass=10**($age_min_mass);
	
    }
    my $name=$unc_file.", ";
    my $scale="1";
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);
    my $pdl_model_spec_min=pdl(@model_spec_min);
    my $pdl_res=pdl(@res_spec);

#    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc);

    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,pdl(@wave_unc),pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);
    }


#    my $end_cat=$ini_cat-1;
#    my $coeffs_cat_sec=$coeffs_cat->slice(":,0:$end_cat");

    
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res,$pdl_rat_now;
}



sub fit_ssp_lin_no_zero {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
    my $sigma_inst=$_[17];
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;
#    my $n_unc=$#flux_unc+1;
    my $n_unc=$#flux_unc+1;

    my $pdl_model_spec_min=zeroes($n_unc);

    my $coeffs=zeroes($nf,3);    
    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc*$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc*$n_mc);


    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];

    my $rsigma=$sigma/$dpix_c;

    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	

	my $pdl_flux_c_ini_now=$pdl_flux_c_ini->slice(",($iii)");
	my $out_spec_pdl=shift_convolve_pdl(pdl(@wave_unc),pdl(@wave_c),$pdl_flux_c_ini_now,0,$sigma_inst,$sigma);
	my ($n_c_out)=$out_spec_pdl->dims;

	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }
#$pdl_flux_c_conv=c
#        $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
#	my $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
#	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
#	my ($n_c_out)=$out_spec_pdl->dims;




    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    my $pdl_dust_spec=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    set($pdl_dust_spec,$i,$j,$dust);
	    my $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
	}
    }
#
# We fit
#
    my $pdl_flux_masked=pdl(@flux_masked);
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);

# Just a linear FIT, without restrictions!
#    ($y_model_now, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;
#    print "$pdl_flux_masked\n $pdl_error\n";
   ($y_model_now, $coeffs) = linfit1d($pdl_flux_masked,$pdl_model,{Weights=>$pdl_error});

#
# We remove the models that are negative
#
    my $nf_new=0;
    $nf_neg=0;
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($C>0) {
	    $nf_new++;
	} else {
	    $nf_neg++;
	}
    }

    my @MOD;
    if ($nf_new>0) {
	while ($nf_neg>0) {
	    my $pdl_model_new=zeroes($n_unc,$nf_new);
	    my $nf_i=0;
	    for ($k=0;$k<$nf;$k++) {
		my $C=$coeffs->at($k,0);
		if ($C>0) {
		    my $t=$pdl_model_new->slice(":,$nf_i");
		    $t .= $pdl_model->slice(":,$k");
		    $MOD[$nf_i]=$k;
		    $nf_i++;
		} else {
		    set($coeffs,$k,0,0);
		}
	    }
	    my $coeffs_new;
	    my $yfit;
	    ($yfit, $coeffs_new) = linfit1d($pdl_flux_masked,$pdl_model_new,{Weights=>$pdl_error});
	    $y_model_now=$yfit;
	    my $nf_i=0;
	    $nf_neg=0;
	    $nf_new=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
		if ($C>0) {
		    my $val=$coeffs_new->at($nf_i,0);
		    $nf_i++;
		    if ($val>0) {
			set($coeffs,$k,0,$val);
			$nf_new++;
		    } else {
			set($coeffs,$k,0,0);
			$nf_neg++;
		    }
		}
	    }
		if ($nf_new==0) {
		    $nf_neg=0;
		}
	}
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    }

    } else {
	$nf_new=$nf;
    }

#
#
#






# Weights=>ones(1) }, $opthash ) ;
##############################
# CHISQ VALUE
	    my $chi=0;
	    my $chi2=0;
	    my $NFREE=0;
	    my @out_spec;
	    my @chi_sec;
	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec[$j]=$y_model_now->at($j,0);
		$model_spec[$j]=$out_spec[$j];
		$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
		$model_spec_min[$j]=$model_spec[$j];
		$chi_sec[$j]=0;
		if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		    $chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    $NFREE++;
		}
	    }
	    
	    my $chi_sq=$chi;
	    if ($NFREE>0) {
		$chi_sq=($chi_sq/($NFREE))**0.5;
	    }

    $chi_sq_min_now=$chi_sq;
    $y_model_end=$y_model_now;
    $y_model_no_mask_end=$y_model_no_mask_now;

    $min_coeffs=$coeffs;

    my $chi_sq=$chi_sq_min_now;
	
    my $norm_C=0;
    my $norm_C_mass=0;

    for ($k=0;$k<$nf;$k++) {
	my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	my $C=$coeffs->at($k,0);
	$norm_C=$norm_C+$C;
	$norm_C_mass=$norm_C_mass+$C*$ml[$k];
	set($coeffs_N,$k,0,$norm_C);
	set($coeffs_NM,$k,0,$norm_C_mass);
	$pdl_model_spec_min=$pdl_model_spec_min+$C*$pdl_model_no_mask->slice(":,($k)");
    }

	
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($norm_C>0) {
	    $age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
	    $met_min=$met_min+$C*$met_mod[$k]/$norm_C;
	    $Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
	    $CN=$C/$norm_C;
	}
	my $C_now=$C*$med_norm;
	if ($norm_C_mass>0) {
	    $age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
	    $met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
	    $Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	}
	
	}
    $age_min=10**($age_min);
    $age_min_mass=10**($age_min_mass);
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);

    my $pdl_res=pdl(@flux_unc)-$pdl_model_spec_min;
    my $pdl_wave_unc=pdl(@wave_unc);
    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,$pdl_wave_unc,pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);

    }
    
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res;
}


sub fit_ssp_lin_no_zero_no_cont {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
    my $sigma_inst=$_[17];
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;
    my $n_unc=$#flux_unc+1;

    my $pdl_model_spec_min=zeroes($n_unc);

    my $coeffs=zeroes($nf,3);    
    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc*$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc*$n_mc);


    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];


    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	

	my $pdl_flux_c_ini_now=$pdl_flux_c_ini->slice(",($iii)");
	my $out_spec_pdl=shift_convolve_pdl(pdl(@wave_unc),pdl(@wave_c),$pdl_flux_c_ini_now,0,$sigma_inst,$sigma);
	my ($n_c_out)=$out_spec_pdl->dims;
	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }



    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    my $pdl_dust_spec=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    set($pdl_dust_spec,$i,$j,$dust);
	    my $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
	}
    }
#
# We fit
#
    my $pdl_flux_masked=pdl(@flux_masked);
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);

    
# Just a linear FIT, without restrictions!
#    ($y_model_now, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;
   ($y_model_now, $coeffs) = linfit1d($pdl_flux_masked,$pdl_model,{Weights=>$pdl_error});

    
    $sigma_mean=sqrt($sigma_inst**2+(5000*($sigma/300000))**2);
    my $smooth_ratio=smooth_ratio(pdl(@flux_unc),$y_model_now,int($sigma_mean));
    $y_model_now=$y_model_now*$smooth_ratio;


#
# We remove the models that are negative
#
    my $nf_new=0;
    $nf_neg=0;
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($C>0) {
	    $nf_new++;
	} else {
	    $nf_neg++;
	}
    }

    my @MOD;
    if ($nf_new>0) {
	while ($nf_neg>0) {
	    my $pdl_model_new=zeroes($n_unc,$nf_new);
	    my $nf_i=0;
	    for ($k=0;$k<$nf;$k++) {
		my $C=$coeffs->at($k,0);
		if ($C>0) {
		    my $t=$pdl_model_new->slice(":,$nf_i");
		    $t .= $pdl_model->slice(":,$k");
		    $MOD[$nf_i]=$k;
		    $nf_i++;
		} else {
		    set($coeffs,$k,0,0);
		}
	    }
	    my $coeffs_new;
	    my $yfit;
	    ($yfit, $coeffs_new) = linfit1d($pdl_flux_masked,$pdl_model_new,{Weights=>$pdl_error});
	    $sigma_mean=sqrt($sigma_inst**2+(5000*($sigma/300000))**2);
	    my $smooth_ratio=smooth_ratio(pdl(@flux_unc),$yfit,int($sigma_mean));
	    $yfit=$yfit*$smooth_ratio;
	    $y_model_now=$yfit;


	    my $nf_i=0;
	    $nf_neg=0;
	    $nf_new=0;
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
		if ($C>0) {
		    my $val=$coeffs_new->at($nf_i,0);
		    $nf_i++;
		    if ($val>0) {
			set($coeffs,$k,0,$val);
			$nf_new++;
		    } else {
			set($coeffs,$k,0,0);
			$nf_neg++;
		    }
		}
	    }
		if ($nf_new==0) {
		    $nf_neg=0;
		}
	}
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    }

    } else {
	$nf_new=$nf;
    }

#
#
#






# Weights=>ones(1) }, $opthash ) ;
##############################
# CHISQ VALUE
	    my $chi=0;
	    my $chi2=0;
	    my $NFREE=0;
	    my @out_spec;
	    my @chi_sec;
	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec[$j]=$y_model_now->at($j,0);
		$model_spec[$j]=$out_spec[$j];
		$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
		$model_spec_min[$j]=$model_spec[$j];
		$chi_sec[$j]=0;		
		if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		    $chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    $NFREE++;
		}
	    }
	    
	    my $chi_sq=$chi;
	    if ($NFREE>0) {
		$chi_sq=($chi_sq/($NFREE))**0.5;
	    }

    $chi_sq_min_now=$chi_sq;
    $y_model_end=$y_model_now;
    $y_model_no_mask_end=$y_model_no_mask_now;

    $min_coeffs=$coeffs;

    my $chi_sq=$chi_sq_min_now;
	
    my $norm_C=0;
    my $norm_C_mass=0;

    for ($k=0;$k<$nf;$k++) {
	my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	my $C=$coeffs->at($k,0);
	$norm_C=$norm_C+$C;
	$norm_C_mass=$norm_C_mass+$C*$ml[$k];
	set($coeffs_N,$k,0,$norm_C);
	set($coeffs_NM,$k,0,$norm_C_mass);
	$pdl_model_spec_min=$pdl_model_spec_min+$C*$pdl_model_no_mask->slice(":,($k)");
    }
    $sigma_mean=sqrt($sigma_inst**2+(5000*($sigma/300000))**2);
    my $smooth_ratio=smooth_ratio(pdl(@flux_unc),$pdl_model_spec_min,int($sigma_mean));
    $pdl_model_spec_min=$pdl_model_spec_min*$smooth_ratio;

	
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($norm_C>0) {
	    $age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
	    $met_min=$met_min+$C*$met_mod[$k]/$norm_C;
	    $Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
	    $CN=$C/$norm_C;
	}
	my $C_now=$C*$med_norm;
	if ($norm_C_mass>0) {
	    $age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
	    $met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
	    $Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	}
	
	}
    $age_min=10**($age_min);
    $age_min_mass=10**($age_min_mass);
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);

    my $pdl_res=pdl(@flux_unc)-$pdl_model_spec_min;
    my $pdl_wave_unc=pdl(@wave_unc);
    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,$pdl_wave_unc,pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);

    }
    
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res;
}


sub fit_ssp_INPUT {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
    my $sigma_inst=$_[17];
    my $coeffs=$_[18];    
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;
    my $n_unc=$#flux_unc+1;

    my $pdl_model_spec_min=zeroes($n_unc);

#    my $coeffs=zeroes($nf,3);    

    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc*$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc*$n_mc);


    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];


    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	

	my $pdl_flux_c_ini_now=$pdl_flux_c_ini->slice(",($iii)");
	my $out_spec_pdl=shift_convolve_pdl(pdl(@wave_unc),pdl(@wave_c),$pdl_flux_c_ini_now,0,$sigma_inst,$sigma);
	my ($n_c_out)=$out_spec_pdl->dims;
	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }



    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    my $pdl_dust_spec=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    set($pdl_dust_spec,$i,$j,$dust);
	    my $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
	}
    }
#
# We fit
#
    my $pdl_flux_masked=pdl(@flux_masked);
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);

# Just a linear FIT, without restrictions!
#   ($y_model_now, $coeffs) = linfit1d($pdl_flux_masked,$pdl_model,{Weights=>$pdl_error});
    my $y_model_now=zeroes($n_unc);
#    print "N_C=$n_c\n";
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	$y_model_now=$y_model_now+$C*$pdl_model->slice(":,($k)");
    }
    
#    $sigma_mean=sqrt($sigma_inst**2+(5000*($sigma/300000))**2);
#    my $smooth_ratio=smooth_ratio(pdl(@flux_unc),$y_model_now,int($sigma_mean));
#    $y_model_now=$y_model_now*$smooth_ratio;
# CHISQ VALUE
    my $chi=0;
    my $chi2=0;
    my $NFREE=0;
    my @out_spec;
    my @chi_sec;
    for ($j=0;$j<$n_unc;$j++) {
	$out_spec[$j]=$y_model_now->at($j,0);
	$model_spec[$j]=$out_spec[$j];
	$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
	$model_spec_min[$j]=$model_spec[$j];
	$chi_sec[$j]=0;		
	if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
	    $chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
	    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
	    $NFREE++;
	}
    }
    
    my $chi_sq=$chi;
    if ($NFREE>0) {
	$chi_sq=($chi_sq/($NFREE))**0.5;
    }

    $chi_sq_min_now=$chi_sq;
    $y_model_end=$y_model_now;
    $y_model_no_mask_end=$y_model_no_mask_now;

    $min_coeffs=$coeffs;

    my $chi_sq=$chi_sq_min_now;
	
    my $norm_C=0;
    my $norm_C_mass=0;

    for ($k=0;$k<$nf;$k++) {
	my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	my $C=$coeffs->at($k,0);
	$norm_C=$norm_C+$C;
	$norm_C_mass=$norm_C_mass+$C*$ml[$k];
	set($coeffs_N,$k,0,$norm_C);
	set($coeffs_NM,$k,0,$norm_C_mass);
	$pdl_model_spec_min=$pdl_model_spec_min+$C*$pdl_model_no_mask->slice(":,($k)");
    }
    $sigma_mean=sqrt($sigma_inst**2+(5000*($sigma/300000))**2);
#    my $smooth_ratio=smooth_ratio(pdl(@flux_unc),$pdl_model_spec_min,int($sigma_mean));
#    $pdl_model_spec_min=$pdl_model_spec_min*$smooth_ratio;

	
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($norm_C>0) {
	    $age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
	    $met_min=$met_min+$C*$met_mod[$k]/$norm_C;
	    $Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
	    $CN=$C/$norm_C;
	}
	my $C_now=$C*$med_norm;
	if ($norm_C_mass>0) {
	    $age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
	    $met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
	    $Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	}
	
	}
    $age_min=10**($age_min);
    $age_min_mass=10**($age_min_mass);
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);

    my $pdl_res=pdl(@flux_unc)-$pdl_model_spec_min;
    my $pdl_wave_unc=pdl(@wave_unc);
    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,$pdl_wave_unc,pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);

    }
    
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res;
}


sub fit_ssp_lin {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $crval=$_[3];
    my $cdelt=$_[4];
    my $crpix=$_[5];
    my $nf=$_[6];
    my $n_c=$_[7];
    my $pdl_flux_c_ini=$_[8];
    my $junk=$_[9];
    my @wave_unc=@$junk;
    my $junk=$_[10];
    my @masked=@$junk;
    my $junk=$_[11];
    my @e_flux_unc=@$junk;
    my $junk=$_[12];
    my @flux_unc=@$junk;
    my @flux_masked;
    my $n_mc=$_[13];
    my $chi_sq_min_now=$_[14];
    my $MIN_CHISQ=$_[15];
    my $plot=$_[16];
    my $k,$j,$i,$iter;
    my $iii,$jj;
    my $iter_max=5;
    my $last_chi=1e12;
#    my $n_unc=$#flux_unc+1;
    my $n_unc=$#flux_unc;

    my $coeffs=zeroes($nf,3);    
    my $coeffs_N=zeroes($nf,1);    
    my $coeffs_NM=zeroes($nf,1);    

    my $coeffs_cat=zeroes($nf+1,$n_mc*$n_mc);

    my $pdl_model_spec_cat=zeroes($n_unc,$n_mc*$n_mc);


    for ($i=0;$i<$n_unc;$i++) {
	$flux_masked[$i]=$flux_unc[$i]*$masked[$i];
    }

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }

    my @wave_c;
    my @dpix_c_val;
    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }

    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];

    my $rsigma=$sigma/$dpix_c;

    my @name;
    my @age_mod,@met_mod;
    my @ml;
    my @flux_c;
    my $pdl_flux_c_conv;


    my @model;
    my @model_no_mask;
    
    my $y_model_end;
    my $y_model_no_mask_end;


    my @MED_NORM;

	my $age_min;
	my $met_min;
	my $Av_min;
	my $age_min_mass;
	my $met_min_mass;
	my $Av_min_mass;

    for ($iii=0;$iii<$nf;$iii++) {
	my $header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini->hdr->{$header};
	my $name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;    
	$name_min =~ s/.dat//;
	my $AGE,$MET;
	my $age,$met;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;	
	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
	$header="NORM".$iii;	
	my $val_ml=$pdl_flux_c_ini->hdr->{$header};
	if ($val_ml!=0) {
	    $ml[$iii]=1/$val_ml;
	} else {
	    $ml[$iii]=1;
	}
	
	my $box=int(3*$rsigma);
	if ($box<3) {
	    $box=3;
	}
	my $kernel=zeroes(2*$box+1);
	my $norm=0;
	$flux_c[$i]=0;
	for ($j=0;$j<2*$box+1;$j++) {
	    $gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	    set($kernel,$j,$gaus);
	    $norm=$norm+$gaus;
	}
	$kernel=$kernel/$norm;
        $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
	my $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	my ($n_c_out)=$out_spec_pdl->dims;
	my @error;

	for ($i=0;$i<$n_unc;$i++) {
	    my $val=$out_spec_pdl->at($i);
	    if ($val eq "nan") {
		$val=0;
	    }
	    $model[$i][$iii]=$val*$masked[$i];
	    $model_no_mask[$i][$iii]=$val;#*$masked[$i];
	    if ($masked[$i]>0) {
		$error[$i]=0.01*abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }

	}
    }



    my $pdl_model=zeroes($n_unc,$nf);
    my $pdl_model_no_mask=zeroes($n_unc,$nf);
    my $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    my $wave_res=$wave_unc[$i]/(1+$redshift);
	    my $dust_rat=A_l(3.1,$wave_res);
	    my $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    my $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    my $val_no_mask=$model_no_mask[$i][$j]*$dust;
	    set($pdl_model_no_mask,$i,$j,$val_no_mask);
	    my $e_val=$error[$i];
	    my $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
	    set($pdl_error,$i,1/(abs($val_now)**2));
	}
    }
#
# We fit
#
    my $pdl_flux_masked=pdl(@flux_masked);
    my $j_mc;
    my $ini_cat=0;
    my $pdl_C_flast=ones($nf);

# Just a linear FIT, without restrictions!
#    ($y_model_now, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;
   ($y_model_now, $coeffs) = linfit1d($pdl_flux_masked,$pdl_model,{Weights=>$pdl_error});
# Weights=>ones(1) }, $opthash ) ;
##############################
# CHISQ VALUE
	    my $chi=0;
	    my $chi2=0;
	    my $NFREE=0;
	    my @out_spec;
	    my @chi_sec;
	    for ($j=0;$j<$n_unc;$j++) {
		$out_spec[$j]=$y_model_now->at($j,0);
		$model_spec[$j]=$out_spec[$j];
		$res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
		$model_spec_min[$j]=$model_spec[$j];
		$chi_sec[$j]=0;
		if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		    $chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		    $NFREE++;
		}
	    }
	    
	    my $chi_sq=$chi;
	    if ($NFREE>0) {
		$chi_sq=($chi_sq/($NFREE))**0.5;
	    }

    $chi_sq_min_now=$chi_sq;
    $y_model_end=$y_model_now;
    $y_model_no_mask_end=$y_model_no_mask_now;

    $min_coeffs=$coeffs;

    my $chi_sq=$chi_sq_min_now;
	
    my $norm_C=0;
    my $norm_C_mass=0;
    for ($k=0;$k<$nf;$k++) {
	my $dust=10**(-0.4*$Av[$k]*$dust_rat);
	my $C=$coeffs->at($k,0);
	$norm_C=$norm_C+$C;
	    $norm_C_mass=$norm_C_mass+$C*$ml[$k];
	set($coeffs_N,$k,0,$norm_C);
	set($coeffs_NM,$k,0,$norm_C);
    }
    

	
    for ($k=0;$k<$nf;$k++) {
	my $C=$coeffs->at($k,0);
	if ($norm_C>0) {
	    $age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
	    $met_min=$met_min+$C*$met_mod[$k]/$norm_C;
	    $Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
	    $CN=$C/$norm_C;
	}
	my $C_now=$C*$med_norm;
	if ($norm_C_mass>0) {
	    $age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
	    $met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
	    $Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	}
	
	}
    $age_min=10**($age_min);
    $age_min_mass=10**($age_min_mass);
    
    my $pdl_age_mod=pdl(@age_mod);
    my $pdl_met_mod=pdl(@met_mod);
    my $pdl_ml=pdl(@ml);
    my $pdl_Av=pdl(@Av);
    my $pdl_model_spec_min=pdl(@model_spec_min);
    my $pdl_res=pdl(@res_spec);

    my $out_ps_now="junk";
    my $title="X=$chi_sq Av=$Av[0] z=$redshift sigma=$sigma";
    if ($plot==1) {
	plot_results($plot,pdl(@wave_unc),pdl($pdl_flux_masked,$pdl_model_spec_min,$pdl_res),$out_ps_now,$title);
    }

    
    return $chi_sq,$pdl_age_mod,$pdl_met_mod,$pdl_ml,$pdl_Av,$coeffs,$coeffs_N,$coeffs_NM,$pdl_model_spec_min,$pdl_res;
}


sub FIT_SINGLE {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;

    my $i,$j,$k;

#
# We renormalize the models
#


    for ($j=0;$j<$n_c_new;$j++) {
	$wave_c_new[$j]=($crval_new+$cdelt_new*($j+1-$crpix_new))*(1+$redshift);
    }
    $dpix_c=$wave_c_new[1]-$wave_c_new[0];
    $rsigma=$sigma/$dpix_c;


    my $CHI_MIN=1e12;

    $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    $kernel=zeroes(2*$box+1);
    $norm=0;
    $flux_c[$i]=0;
    for ($j=0;$j<2*$box+1;$j++) {
	$gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
    }
    $kernel=$kernel/$norm;
    
    $pdl_flux_c_conv = conv2d $pdl_flux_c_ini_new,$kernel;
    

    for ($iii=0;$iii<$nf_new;$iii++) {
	$header="NAME".$iii;
	$name[$iii]=$pdl_flux_c_ini_new->hdr->{$header};
	$name_min=$name[$iii];
	$name_min =~ s/spec_ssp_//;
	$name_min =~ s/.spec//;
	($AGE,$MET)=split("_",$name_min);
	if ($AGE =~ "Myr") {
	    $age=$AGE;
	    $age =~ s/Myr//;
	    $age=$age/1000;
	} else {
	    $age=$AGE;
	    $age =~ s/Gyr//;
	}
	$met=$MET;
	$met =~ s/z/0\./;

	$age_mod[$iii]=$age;
	$met_mod[$iii]=$met;
#****** CONVOLUTION
	

	$pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");

	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c_new), $pdl_flux_c);
	($n_c_out)=$out_spec_pdl->dims;
	

	
	
	for ($i=0;$i<$n_unc;$i++) {
	    if ($masked[$i]>0) {
		$error[$i]=0.01/abs($e_flux_unc[$i]);
	    } else {
		$error[$i]=0;
	    }
	    $val=$out_spec_pdl->at($i);
	    $model[$i][$iii]=$val;#*$masked[$i];
#	    if ($NITER==0) {
		for ($j=0;$j<$nline;$j++) {
		    if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
#		    $error[$i]=(30/(1+$INTER))*$e_flux_unc[$i];
		    $error[$i]=0.1/abs($e_flux_unc[$i]);
#			$error[$i]=10*$e_flux_unc[$i];
		    }
		}
	    }
#	}
	
    }



    $pdl_model=zeroes($n_unc,$nf_new);
    $pdl_error=zeroes($n_unc,$nf_new);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf_new;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
#	    if ($j==0) {
		$wave_res=$wave_unc[$i]/(1+$redshift);
		$dust_rat=A_l(3.1,$wave_res);
		$dust=10**(-0.4*$Av[$j]*$dust_rat);  
		$a_dust[$i]=$dust;
#	    }
	    $val=$model[$i][$j]*$a_dust[$i];
	    set($pdl_model,$i,$j,$val);
	    $e_val=$error[$i];
	    set($pdl_error,$i,$j,1);
	}
    }
    $pdl_flux_masked=pdl(@flux_masked);


#    ($yfit, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;
    #
    # We normalize the models
    #
 

#    my $nx1=int(($wave_norm-$w_wave_norm-$crval3)/$cdelt3);
#    my $nx2=int(($wave_norm+$w_wave_norm-$crval3)/$cdelt3);
#    my $pdl_sec_cen=$pdl_flux_masked->slice("$nx1:$nx2");
#    my $SUM=sum($pdl_sec_cen);
    
    my $pdl_chi=zeroes($nf_new);
    my $pdl_e=pdl(@error);
    $J_MIN=-1;
   for ($j=0;$j<$nf_new;$j++) {
       my $pdl_model_tmp=$pdl_model->slice(":,$j");
       ($yfit, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model_tmp,$pdl_e;


       my $pdl_chi_now=$pdl_masked*(($pdl_flux_masked-$yfit)**2)*$pdl_e;
       my $chi_sq=sum($pdl_chi_now);
       $chi_sq=($chi_sq/($n_unc-1))**0.5;
       set($pdl_chi,$j,$chi_sq);

       if ($CHI_MIN>$chi_sq) {
	   $CHI_MIN=$chi_sq;
	   $J_MIN=$j;
	   $NORM_J=$NORM;


#	   @out_spec=list($pdl_model_tmp);
#	   @model_spec_min=list($pdl_model_tmp);
	   @out_spec=list($yfit);
	   @model_spec_min=list($yfit);
	   my $pdl_res=pdl(@flux_unc)-$pdl_model_tmp;
	   @res_spec=list($pdl_res);
	   $JJ=0;
	   for ($jj=0;$jj<$n_unc;$jj++) {
	       if ($masked[$jj]>0) {
		   $res_clean[$JJ]=$res_spec[$jj];
		   $flux_clean[$JJ]=$flux_unc[$jj];
		   $JJ++;
	       }
	   }




	if ($plot==1) {

	    pgbegin(0,$dev_plot_single,1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.2);           # Set character height
	    pgenv($wave_unc[0],$wave_unc[$n_unc-1],$y_min,$y_max,0,0);
	    pgsch(1.2);           # Set character height
	    pglabel("Wavelength","Counts","SSP X=$CHI_MIN T=$age_mod[$J_MIN] Z=$met_mod[$J_MIN] z=$redshift sigma=$sigma Av=@Av");
	    pgsci(1);
	    pgline($n_unc,\@wave_unc,\@flux_unc);    
	    pgsci(2);
	    pgline($n_unc,\@wave_unc,\@out_spec);    
	    pgsci(3);
	    pgline($n_unc,\@wave_unc,\@flux_masked);    
	    pgsci(1);
	    pgclose;
	    pgend;
	}

       }

   }
#    <stdin>;
#    $pdl_model_tmp=$pdl_model->slice(":,$J_MIN");
#    $pdl_model_tmp=$pdl_model_tmp*$NORM_J;



#    print "SSP CHI=$CHI_MIN $age_mod[$J_MIN] $met_mod[$J_MIN]\n";

    if ($CHI_MIN<$MIN_CHISQ) {    
	$MIN_CHISQ=$CHI_MIN;
	@Av_MIN=@Av;
	$age_min=$age_mod[$J_MIN];
	$met_min=$met_mod[$J_MIN];

    }
    return($CHI_MIN);
}



################ FIT_NEW #########################


sub FIT_Z_SIGMA {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
    my $k,$j,$i,$iter;

    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
	}
	
    }
    $dpix_c_val[0]=$dpix_c_val[1];
    $dpix_c=$wave_c[1]-$wave_c[0];
    $rsigma=$sigma/$dpix_c;

# CONVOLUTION
    $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    $kernel=zeroes(2*$box+1);
    $norm=0;
    for ($j=0;$j<2*$box+1;$j++) {
	$gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
    }
    $kernel=$kernel/$norm;

    $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;
# END CONVOLUTION

    $pdl_dust=zeroes($n_unc,$nf);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    $wave_res=$wave_unc[$i]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    if ($Av[$j]!=0) {
		$dust=10**(-0.4*$Av[$j]*$dust_rat);  
		set($pdl_dust,$i,$j,$dust);
	    } else {
		set($pdl_dust,$i,$j,1);
	    }
	}
    }


    my $model_now=zeroes($n_unc);    
    for ($iii=0;$iii<$nf;$iii++) {    
	$pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");
	my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
	$C=$coeffs->at($iii,0);
	my $dust_spec=$out_spec_pdl*$pdl_dust->slice(",$iii");
	$model_now=$model_now+$C*$dust_spec;
    }

    $chi=0;
    $chi2=0;
    $NFREE=0;
    for ($j=0;$j<$n_unc;$j++) {
	if (($wave_unc[$j]>$MIN_W)&&($wave_unc[$j]<$MAX_W)) {
	    $out_spec[$j]=$model_now->at($j,0);
	    if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)) {
		$chi_j[$j]=$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/abs(1.5*$out_spec[$j]);
		$chi=$chi+$chi_j[$j];
		$NFREE++;
	    }
	}
    }

   $chi_sq=$chi;
    if ($NFREE>0) {
	$chi_sq=($chi_sq/($NFREE))**0.5;
    }
    return $chi_sq;
}





sub FIT_simple {
    my $redshift=$_[0];
    my $sigma=$_[1];
    my $Av_NOW=$_[2];
    my @Av=@$Av_NOW;
#    print "AV FIT = '@AV' '$Av_NOW' '$_'\n";
    my $k,$j,$i,$iter;

    for ($i=0;$i<$nf;$i++) {
	if ($Av[$i]<0) {
	    $Av[$i]=0;
	}
    }


    for ($j=0;$j<$n_c;$j++) {
	$wave_c[$j]=($crval+$cdelt*($j+1-$crpix))*(1+$redshift);
	if ($j>0) {
	    $dpix_c_val[$j]=$wave_c[$j]-$wave_c[$j-1];	
#	    print "$j $dpix_c_val[$j] $crval $cdelt\n";
	}
	
    }
    $dpix_c_val[0]=$dpix_c_val[1];
$dpix_c=$wave_c[1]-$wave_c[0];

#$dpix_c=median($wave_c[1]-$wave_c[0];

#($dpix_c,$dpix_max)=minmax(@dpix_c_val);
#$dpix_c=median(@dpix_c_val);
#print "DPIX= $dpix_c\n"; exit;

$rsigma=$sigma/$dpix_c;

#    print "$sigma $dpix_c\n";

for ($iii=0;$iii<$nf;$iii++) {
    $header="NAME".$iii;
    $name[$iii]=$pdl_flux_c_ini->hdr->{$header};
    $name_min=$name[$iii];
    $name_min =~ s/spec_ssp_//;
    $name_min =~ s/.spec//;    
    $name_min =~ s/.dat//;
    ($AGE,$MET)=split("_",$name_min);
    if ($AGE =~ "Myr") {
	$age=$AGE;
	$age =~ s/Myr//;
	$age=$age/1000;
    } else {
	$age=$AGE;
	$age =~ s/Gyr//;
    }
    $met=$MET;
    $met =~ s/z/0\./;

    $age_mod[$iii]=$age;
    $met_mod[$iii]=$met;
    $header="NORM".$iii;
#    $ml[$iii]=$pdl_flux_c_ini->hdr->{$header};

 $val_ml=$pdl_flux_c_ini->hdr->{$header};
    if ($val_ml!=0) {
	$ml[$iii]=1/$val_ml;
    } else {
	$ml[$iii]=1;
    }




#    print "$name[$iii] $age_mod[$iii] $met_mod[$iii]\n";
    
    
    $box=int(3*$rsigma);
    if ($box<3) {
	$box=3;
    }
    $kernel=zeroes(2*$box+1);
    $norm=0;
    $flux_c[$i]=0;
    for ($j=0;$j<2*$box+1;$j++) {
	$gaus=exp(-0.5*((($j-$box)/$rsigma)**2));    
	set($kernel,$j,$gaus);
	$norm=$norm+$gaus;
    }
    $kernel=$kernel/$norm;
    

    $pdl_flux_c_conv = conv2d $pdl_flux_c_ini,$kernel;



    $pdl_flux_c = $pdl_flux_c_conv->slice(",$iii");


     my $out_spec_pdl = interpol(pdl(@wave_unc), pdl(@wave_c), $pdl_flux_c);
    ($n_c_out)=$out_spec_pdl->dims;
#
# We store the result!
#
    for ($i=0;$i<$n_unc;$i++) {
#	print "$i/$n_c\n";
	$val=$out_spec_pdl->at($i);
	if ($val eq "nan") {
	    $val=0;
	}
	$model[$i][$iii]=$val*$masked[$i];
	$model_no_mask[$i][$iii]=$val;#*$masked[$i];

	if ($masked[$i]>0) {
	    $error[$i]=0.01*abs($e_flux_unc[$i]);
#	    $error[$i]=abs($e_flux_unc[$i]);
	} else {
	    $error[$i]=0;
	}


	
#	if ($NITER==0) {
	    for ($j=0;$j<$nline;$j++) {
		if (abs($wave_unc[$i]-$w_eline[$j]*(1+$redshift))<(3*$rsigma)) {
#		    $error[$i]=0.1*$e_flux_unc[$i];
#		    $error[$i]=30;#*$e_flux_unc[$i];
		}
	    }
#	}

#		    print "$i $error[$i]\n";

#	print "$i $wave_unc[$i] $error[$i] $e_flux_unc[$i]\n";

    }
 #   print "$iii/$nf files\n";
}

    $pdl_model=zeroes($n_unc,$nf);
#    $pdl_error=zeroes($n_unc,$nf);
    $pdl_error=zeroes($n_unc);
    my $pdl_masked=pdl(@masked);
    for ($j=0;$j<$nf;$j++) {
	for ($i=0;$i<$n_unc;$i++) {
	    $wave_res=$wave_unc[$i]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $dust=10**(-0.4*$Av[$j]*$dust_rat);  
	    $val=$model[$i][$j]*$dust;
	    set($pdl_model,$i,$j,$val);
	    $e_val=$error[$i];
	    $val_now=$e_flux_unc[$i];
	    if ($val_now==0) {
		$val_now=1;
	    }
#	    set($pdl_error,$i,$j,1/(abs($val_now)**2));
	    set($pdl_error,$i,1/(abs($val_now)**2));
	}
    }
#
# We fit
#
#    print "We fit with Av=$Av ";
    $pdl_flux_masked=pdl(@flux_masked);
    ($yfit, $coeffs) = my_linfit1d $pdl_flux_masked,$pdl_model,$pdl_error;


#
# We remove the models that are negative
#
    my $nf_new=0;
    $nf_neg=0;
    for ($k=0;$k<$nf;$k++) {
	$C=$coeffs->at($k,0);
	if ($C>0) {
	    $nf_new++;
	} else {
	    $nf_neg++;
	}
    }

#    print "\n";

   
    if ($nf_new>0) {
	while ($nf_neg>0) {
#	    print "PASO $nf_new $nf_neg \n";
#	    if ($nf_new>0) {
		my $pdl_model_new=zeroes($n_unc,$nf_new);
#		my $pdl_error_new=zeroes($n_unc,$nf_new);
		my $pdl_error_new=zeroes($n_unc);
#	    }
	    $nf_i=0;
#	    print "C= ";
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
#		print " $C ";
#		print "$C $k PASO\n";
		if ($C>0) {
		    for ($i=0;$i<$n_unc;$i++) {
			$val=$pdl_model->at($i,$k);
			set($pdl_model_new,$i,$nf_i,$val);
#			print "$i PASO\n";
#			$e_val=$pdl_error->at($i,$k);
			$e_val=$pdl_error->at($i);
#			set($pdl_error_new,$i,$nf_i,$e_val);		    
			set($pdl_error_new,$i,$e_val);		    
		    }
		    $MOD[$nf_i]=$k;
		    $nf_i++;
		} else {
		    set($coeffs,$k,0,0);
		}
	    }
#	    print "\n";
#	    print "PASO\n";
	    ($yfit, $coeffs_new) = my_linfit1d $pdl_flux_masked,$pdl_model_new,$pdl_error_new;	
#	    print "FIT $coeffs_new\n";
	    
	    $nf_i=0;
	    $nf_neg=0;
	    $nf_new=0;
#	    print "C= ";
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
#		print "$C ";
		if ($C>0) {
		    $val=$coeffs_new->at($nf_i,0);
		    $nf_i++;
		    if ($val>0) {
			set($coeffs,$k,0,$val);
			$nf_new++;
		    } else {
			set($coeffs,$k,0,0);
			$nf_neg++;
		    }
#		    print "$val , ";
		}
	    }
#	    print "\n";
		if ($nf_new==0) {
		    $nf_neg=0;
		}
	}
#	    print "$nf_new $nf_neg ";
#	    print "C= ";
	    for ($k=0;$k<$nf;$k++) {
		$C=$coeffs->at($k,0);
#		print " $C ";
	    }
#	print "\n";

    } else {
	$nf_new=$nf;
	$pdl_model_new=$pdl_model;
	$pdl_error_new=$pdl_error;
    }



#    print "PASO\n";

#
# We determine the chi**2
#
    @J=$yfit->dims;
    @K=$coeffs->dims;
#    print "DONE ( $coeffs)\n"; <stdin>;
    $chi=0;
    $chi2=0;
    $NFREE=0;
    for ($j=0;$j<$n_unc;$j++) {
	$out_spec[$j]=($yfit->at($j,0));#
	$chi_sec[$j]=0;
	    if (($flux_unc[$j]!=0)&&($out_spec[$j]!=0)&&($e_flux_unc[$j]!=0)) {
		#print "$j $flux_unc[$j]-$out_spec[$j]\n";
#0000000000000000000000000000000
#		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs(1.5*$out_spec[$j]);
		$chi=$chi+$masked[$j]*(($flux_masked[$j]-$out_spec[$j])**2)/abs(1.5*$out_spec[$j]);
#		$chi2=$chi2+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($e_flux_unc[$j]);
		if ($have_error==0) {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/abs($out_spec[$j]);
		} else {
		    $chi_sec[$j]=$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)/($e_flux_unc[$j])**2;
		}
#		$chi=$chi+$masked[$j]*(($flux_unc[$j]-$out_spec[$j])**2)*$error[$j];
		$NFREE++;
	    }
    }

    $chi_sq=$chi;
    if ($NFREE>0) {
	$chi_sq=($chi_sq/($NFREE))**0.5;
    }
#    $chi2=($chi2/($NFREE))**0.5;

#    $chi_sq=(($chi/($n_unc-$nf))**0.5);
#    $chi_sq=(($chi/($n_unc-1))**0.5);
    if ($chi_sq<$MIN_CHISQ) {
#    print "PASO $chi_sq $MIN_CHISQ\n";
	$MIN_CHISQ=$chi_sq;
	@Av_MIN=@Av;
	for ($j=0;$j<$n_unc;$j++) {
	    $model_spec[$j]=0;
	    $wave_res=$wave_unc[$j]/(1+$redshift);
	    $dust_rat=A_l(3.1,$wave_res);
	    $norm_C=0;
	    $norm_C_mass=0;
	    for ($k=0;$k<$nf;$k++) {
		$dust=10**(-0.4*$Av[$k]*$dust_rat);  
		$C=$coeffs->at($k,0);
		$norm_C=$norm_C+$C;
		$norm_C_mass=$norm_C_mass+$C*$ml[$k];
		$model_spec[$j]=$model_spec[$j]+$C*$model_no_mask[$j][$k]*$dust;
		$model_spec_min[$j]=$model_spec[$j];
	    }
	    $res_spec[$j]=$flux_unc[$j]-$model_spec[$j];	    
#	    print "$wave_unc[$j] $flux_unc[$j] $model_spec[$j] $res_spec[$j]\n";
	}
	$age_min=0;
	$met_min=0;
	$Av_min=0;
	$age_min_mass=0;
	$met_min_mass=0;
	$Av_min_mass=0;
	open(C,">coeffs.out");
	print C "ID   AGE     MET    COEFF   Norm.Coeff  M/L   AV\n";
	print "---------------------------------------------\n";
	print "ID AGE     MET    COEFF   Norm.Coeff   M/L    AV\n";
	print "---------------------------------------------\n";
	for ($k=0;$k<$nf;$k++) {
	    $C=$coeffs->at($k,0);
	    if ($norm_C>0) {
		$age_min=$age_min+$C*log10($age_mod[$k])/$norm_C;
		$met_min=$met_min+$C*$met_mod[$k]/$norm_C;
		$Av_min=$Av_min+$C*$Av[$k]/$norm_C;	    
		$CN=$C/$norm_C;
	    }
#	    if ($C==0) {
#		$d_Av[$k]=0;
#	    }
#	    print C "$k $age_mod[$k] $met_mod[$k] $C $CN $ml[$k] $Av[$k]\n";
	    printf(C "%-2d %-7.4f %-7.4f %-7.4f %-7.4f %-10.2f %-4.2f\n",$k,$age_mod[$k],$met_mod[$k],$C,$CN,$ml[$k],$Av[$k]);
	    if ($C>0) {
		printf("%-2d %-7.4f %-7.4f %-7.4f %-7.4f %-10.2f %-4.2f\n",$k,$age_mod[$k],$met_mod[$k],$C,$CN,$ml[$k],$Av[$k]);
	    }
	    if ($norm_C_mass>0) {
		$age_min_mass=$age_min_mass+$C*log10($ml[$k]*$age_mod[$k])/$norm_C_mass;
		$met_min_mass=$met_min_mass+$C*$ml[$k]*$met_mod[$k]/$norm_C_mass;
		$Av_min_mass=$Av_min_mass+$C*$ml[$k]*$Av[$k]/$norm_C_mass;
	    }
#		print "COEFF=$k $C $norm_C_mass $age_min $age_mod[$k] $met_min $met_mod[$k] $ml[$k] $age_min_mass $met_min_mass\n";

	}
	print "---------------------------------------------\n";
	close(C);
	$age_min=10**($age_min);
	$age_min_mass=10**($age_min_mass);
	
    }
    $name=$unc_file.", ";
    $scale="1";
#    print "*****PASO 1***** PLOT = $plot, N_UNC = $n_unc, $nc\n";
    
    if ($plot==1) {

	pgbegin(0,$dev_plot,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($min_wave,$max_wave,$y_min,$y_max,0,0);
	pgsch(0.9);           # Set character height
#	pglabel("Wavelength","Counts","X=$chi_sq ($chi_joint) T=$age_min ($age_min_mass) Z=$met_min ($met_min_mass) Av=@Av z=$redshift sigma=$sigma");
	my $mean_Av=mean(@Av);
	pglabel("Wavelength","Counts","X=$chi_sq ($chi_joint) T=$age_min ($age_min_mass) Z=$met_min ($met_min_mass) Av=$mean_Av z=$redshift sigma=$sigma");
#	pglabel("Wavelength","Counts","X=$chi_sq ($chi_joint) T=$age_min ($age_min_mass) Z=$met_min ($met_min_mass) z=$redshift sigma=$sigma");
	pgsch(0.5);
	pgsci(1);
#	print "**** PASO NOW ****\n";


	pgline($n_unc,\@wave_unc,\@flux_unc);    
	pgsci(2);
	pgline($n_unc,\@wave_unc,\@out_spec);    
	pgsci(8);
	pgline($n_unc,\@wave_unc,\@model_spec_min);    
	pgsci(5);
	pgline($n_unc,\@wave_unc,\@res_spec);    

	if ($nc>0) {
	    pgsci(6);
	    pgline($nc,\@wave_clean,\@res_clean);    
	}
	pgsci(1);
	pgsci(13);
	pgsls(2);
	$j_plot=0;



	for ($j=0;$j<$nf;$j++) {
	    $C=$coeffs->at($j,0);
	    if ($C>0) {
		pgsci(13);
		my $model_now=$pdl_model->slice(",($j)");

		$model_now=$C*$model_now;
		@a_model_now=list($model_now);
		pgline($n_unc,\@wave_unc,\@a_model_now);    

	    }
	    if ($norm_C != 0) {
		$CN=$C/$norm_C;
	    } else {
		$CN=$C;
	    }
	    if ($C>0) {
		$pos_x=$min_wave*1.02;
		$pos_y=$y_max-0.05*($j_plot+1)*($y_max-$y_min);
		pgsch(0.8);           # Set character height
		if ($C>0) {
		    pgsci(1);
		} else {
		pgsci(15);
		}
		$CN=apr($CN);
		pgptxt($pos_x,$pos_y,0,0,"$CN $age_mod[$j] $met_mod[$j]");
		$j_plot++;
#	    print "$pos_x $pos_y\n"; <stdin>;
	    }
		pgsch(1.2);           # Set character height

	}
	pgsls(1);


	pgsci(1);


	pgclose;
	pgend;
#    print "Press Enter"; <stdin>;
    }


#
#
#





    return $chi_sq;
}


sub get_time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $time="# TIME ".$sec." ".$min." ".$hour." ".$mday." ".$mon." ".$year." ".$wday." ".$yday." ".$isdst;
    return($time);
}


sub fit_elines_rnd() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_mc=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];
    my @a_out;    
#    my $SCALE_INI=3;
    my $SCALE=$SCALE_INI;
#    print "N.MOD = $n_mod, N.MC=$n_mc\n";

    my @dims=$pdl_flux->dims;
    my $nx,$ny;
    if ($#dims==0) {
	$nx=$dims[0];
	$ny=1;
    }

    if ($#dims==1) {
	$nx=$dims[0];
	$ny=$dims[1];
    }
#    print "$nx,$ny\n";
#    my $pdl_model=zeroes($nx,$ny);
#    my $pdl_masked=ones($nx);

    my $i;
    #
    #

    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0=zeroes($nx);
    my $pdl_model_cont_0=zeroes($nx);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
	if ($type[$i] eq "poly1d") {
	    $pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
	}
    }


#    my $pdl_grandom=grandom($nx);
 #   my $pdl_flux_fit=$pdl_flux+$pdl_grandom*$pdl_e_flux;

#    my $pdl_chi_now=$pdl_masked*(($pdl_flux_fit-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
#    print "CHI_sq [GUESS]=$chi_sq\n";  
	$chi_sq_ini=$chi_sq;
	@a_out=copy_a($n_mod,\@a);




    #############################
    # TEST OUTPUT!
    #
    open(OUT,">test_fit_elines_rnd.spec");
    for ($i=0;$i<$nx;$i++) {
	my $w=$pdl_wave->at($i);
	my $f=$pdl_flux->at($i);
	my $f_mod=$pdl_model_0->at($i)*($pdl_masked->at($i));
	my $f_res=$f-$f_mod;
	my $f_cont=$pdl_model_cont_0->at($i);
	print OUT "$w $f $f_mod $f_red $f_cont\n";
    }
    close(OUT);
    $call="plot_out_fit_mod.pl test_fit_elines_rnd.spec 1/xs 0 2.5e4 6450 6800 &";
    system($call);
#    print "###################################\n";
#    <stdin>;
    #
    # TEST OUTOUT
    ###############################

    my $i,$j,$ii;
    my @a_now;

#    srand(localtime);
	$pdl_rnd=random($n_mod*9*$n_mc);

    for ($ii=0;$ii<$n_mc;$ii++) {
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			$a_now[$i][$j]=$a0[$i][$j]+($a1[$i][$j]-$a0[$i][$j])*$rnd;	
		    } else {
			$k=$link[$i][$j]-1;
			$method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a[$i][$j];
		}	    
	    }
	}

    # MC test
    $pdl_model_0=zeroes($nx);
    $pdl_model_cont_0=zeroes($nx);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
	if ($type[$i] eq "poly1d") {
	    $pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
	}
    }
#    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);

    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
	if ($chi_sq<$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
	    $chi_sq_ini=$chi_sq;


    #############################
    # TEST OUTPUT!
    #
    open(OUT,">test_fit_elines_rnd.spec");
    for ($i=0;$i<$nx;$i++) {
	my $w=$pdl_wave->at($i);
	my $f=$pdl_flux->at($i);
		my $f_mod=$pdl_model_0->at($i)*($pdl_masked->at($i));
	my $f_res=$f-$f_mod;
	my $f_cont=$pdl_model_cont_0->at($i);
	print OUT "$w $f $f_mod $f_red $f_cont\n";
    }
    close(OUT);
#    print "iter= $ii/$n_mc  CHI_sq=$chi_sq / $chi_sq_ini SCALE = $SCALE\n";	    
#	    $call="plot_out_fit_mod.pl test_fit_elines_rnd.spec 1/xs 0 2.5e4 6300 7000 ";
    $call="plot_out_fit_mod.pl test_fit_elines_rnd.spec 1/xs 0 2.5e4 6450 6800 &";
    system($call);
    #
    # TEST OUTOUT
	    ###############################


	    my $iii;
	    for ($iii=0;$iii<$n_mod-1;$iii++) {
		if (($type[$iii] eq "eline") and ($link[$iii][1]==-1)) {
		    $a0[$iii][1]=$a_out[$iii][1]/(1+1/$SCALE);
		    $a1[$iii][1]=$a_out[$iii][1]*(1+$SCALE);
		}
	    }


	    $SCALE=0.75*$SCALE;
	    if ($SCALE>0.001) {
		$ii=0;
	    } else {
		$ii=$n_mc;
	    }








	} else {

	#    $SCALE=1.1*$SCALE;

	}
#else {	    
#	    if (($def==1)&&($chi_sq>2*$chi_sq_ini)) {
#	    if ($def==3) {
#
#		}
#	    }
#	}

	   
    }
    $chi_sq_now=$chi_sq_ini;





    return @a_out;
}




sub get_initial_guess() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $i_now=$_[2];
    my $junk=$_[3];
    my @type=@$junk;
    my $junk=$_[4];
    my @a=@$junk;
    my @dims=$pdl_wave->dims();
    my $nx=$dims[0];
    my $i=0;
#    my $pdl_out=zeroes($nx);
    my $j=0;
    my $guess;

    if ($type[$i_now] eq "eline") {
	$guess=$a[$i_now][1];
	$min_delta_wave=1e12;
	for ($i=0;$i<$nx;$i++) {
	    my $w=$pdl_wave->at($i);	    
	    my $speed_of_light=299792.458;
	    my $factor=(1+$a[$i_now][3]/$speed_of_light+$a[$i_now][5]);
	    my $delta_wave=abs($w-$a[$i_now][0]*$factor);
	    if ($min_delta_wave>$delta_wave) {
		$min_delta_wave=$delta_wave;
		$guess=$pdl_flux->at($i);
		$guess=$guess*($a[$i_now][2]*((2*3.1416)**0.5));
	    }
	}
    }

    if ($type[$i_now] eq "poly1d") {
	$guess=$a[$i_now][0];
	@stats=stats($pdl_flux);
	$guess=$stats[2];
    }
    return $guess;
}



sub fit_elines_grad_rnd_lin() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_mc=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];

    my @a_out=copy_a($n_mod,\@a);

    my $SCALE=$SCALE_INI;
    my $max_time=$_[16];

    my @dims=$pdl_flux->dims;
    my $time1=get_seconds();

    my $nx,$ny;
    if ($#dims==0) {
	$nx=$dims[0];
	$ny=1;
    }

    if ($#dims==1) {
	$nx=$dims[0];
	$ny=$dims[1];
    }

    my $i;

    my @e_stats=stats($pdl_e_flux);
    $pdl_e_flux->inplace->setvaltobad(0);
    $pdl_e_flux->inplace->setbadtoval($e_stats[0]); 

    my $pdl_model=zeroes($nx);
    my $pdl_model_cont=zeroes($nx);
    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0,$coeffs_0;#=zeroes($nx);
    my $pdl_model_cont_0=zeroes($nx);



    my $n_mod_free=0;
    my $i,$j;
    my $pdl_model_n;
    my $cont;
    for ($i=0;$i<$n_mod;$i++) {
	if ($type[$i] eq "eline") {
	    if ($ia[$i][1]==1) {
		if ($link[$i][1]==-1) {
		    my $pdl_tmp=create_single_model_one($pdl_wave,$i,\@type,\@a);
		    if ($n_mod_free==0) {
			$pdl_model_n=$pdl_tmp;
		    } else {
			$pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
		    }
		    $n_mod_free++;
		} else {
		    my $pdl_tmp=create_single_model_one($pdl_wave,$i,\@type,\@a);
		    my $ii=$link[$i][1]-1;
		    my $t=$pdl_model_n->slice(":,$ii");
		    $pdl_tmp=$pdl_tmp*$a0[$i1][1];
		    $t .= $t+$pdl_tmp;
		}
	    }
	}
	if ($type[$i] eq "poly1d") {	    
	    for ($j=0;$j<9;$j++) {
		if ($ia[$i][$j]==1) {
		    my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
		    $pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
		    $n_mod_free++;
		} else {
		    my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
		    $cont=$cont+$a[$i][$j]*$pdl_tmp;
		}
	    }
	}
   }

#    print "$pdl_model_n\n";
#    print "N_MOD_FREE=$n_mod_free, N_MOD=$n_mod\n";

#    my $pdl_model_n=zeroes($nx,$n_mod);
 #   for ($i=0;$i<$n_mod;$i++) {
#	my $pdl_tmp=create_single_model_one($pdl_wave,$i,\@type,\@a);
#	my $t=$pdl_model_n->slice(":,$i");
#	$t .= $pdl_tmp;
 #   }

#    my $pdl_flux_fit=$pdl_flux-$cont;
	my $pdl_grandom=grandom($nx);
	my $pdl_flux_fit=$pdl_flux-$cont+$pdl_grandom*$pdl_e_flux;
    @dim_n=$pdl_model_n->dims();
    if ($dim_n[1]>1) {
	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
    } else {
	my $pdl_rat=$pdl_flux_fit/$pdl_model_n;
	$pdl_rat->inplace->setnantobad();
	$pdl_rat->inplace->setbadtoval(0);
	my @stats=stats($pdl_rat);
	$coeffs_0=ones(1);
	$pdl_model_0=$stats[2]*$pdl_model_n;
	$coeffs_0=$stats[2]*$coeffs_0;
#	print "@stats\n";
    }

#    my $pdl_grandom=grandom($nx);
#    my $pdl_flux_fit=$pdl_flux-$cont+$pdl_grandom*$pdl_e_flux;
#   ($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>1/$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
    $pdl_model_0=$pdl_model_0+$cont;
#    print "$coeffs_0\n";
#    print "$pdl_model_0\n";

    my $n_mod_free=0;
    my $cont;
    for ($i=0;$i<$n_mod;$i++) {
	if ($type[$i] eq "eline") {
	    if ($ia[$i][1]==1) {
		if ($link[$i][1]==-1) {
		    $a[$i][1]=$coeffs_0->at($n_mod_free);
		    $n_mod_free++;
		} else {
		    my $k=$link[$i][1]-1;
		    my $method=$a1[$i][1];
		    if ($method==0) {
			$a[$i][1]=$a[$k][1]+$a0[$i][1];
		    } else {
			$a[$i][1]=$a[$k][1]*$a0[$i][1];
		    }			    
		}
	    }
	}
	if ($type[$i] eq "poly1d") {
	    for ($j=0;$j<9;$j++) {
		if ($ia[$i][$j]==1) {
		    my $val=$coeffs_0->at($n_mod_free);
		    $a[$i][$j]=$val;
		    $n_mod_free++;
		}  else {
		    my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
		    $cont=$cont+$a[$i][$j]*$pdl_tmp;
		}
	    }
	}
    }

    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);

    #$pdl_model=$pdl_model_0;


    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


#    print "";

    my $i,$j,$ii;
    my @a_now;


    # We start the fitting loop!
    $pdl_rnd=grandom($n_mod*9*$n_mc);
    $pdl_rnd_lin=random($n_mod*9*$n_mc);
    my $ii=0;
    my $pdl_flux_in=$pdl_flux;
    while ($ii<$n_mc) {
#	print "$ii/$n_mc\n";

	#
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			if ($j==3) {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$rnd*($a1[$i][$j]-$a0[$i][$j])/5;
			} else {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			}
			if ($a_now[$i][$j]<$a0[$i][$j]) {
			    $a_now[$i][$j]=$a0[$i][$j];
			}
			if ($a_now[$i][$j]>$a1[$i][$j]) {
			    $a_now[$i][$j]=$a1[$i][$j];
			}
		    } else {
			my $k=$link[$i][$j]-1;
			my $method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
#		print "$a_now[$i][$j] $ia[$i][$j] $a0[$i][$j] $a1[$i][$j] $link[$i][$j]\n";
	    }
	}


	
	# MC test
#	$pdl_model_0=zeroes($nx);
#	$pdl_model_cont_0=zeroes($nx);
#	for ($i=0;$i<$n_mod;$i++) {
#	    my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
#	    $pdl_model_0=$pdl_model_0+$pdl_tmp;
#	    if ($type[$i] eq "poly1d") {
#		$pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
#	    }
#	}

	my $n_mod_free=0;
	my $i1,$j1;
	my $pdl_model_n;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			if ($n_mod_free==0) {
			    $pdl_model_n=$pdl_tmp;
			} else {
			    $pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			}
			$n_mod_free++;
		    } else {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			my $ii=$link[$i1][1]-1;
			my $t=$pdl_model_n->slice(":,$ii");
			$pdl_tmp=$pdl_tmp*$a0[$i1][1];		
			$t .= $t +$pdl_tmp;		   
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {	    
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
			$pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			$n_mod_free++;
		    }
		}
	    }
	}


	my $pdl_grandom=grandom($nx);
	my $pdl_flux_fit=$pdl_flux-$cont+$pdl_grandom*$pdl_e_flux;
    @dim_n=$pdl_model_n->dims();
    if ($dim_n[1]>1) {
	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>1/$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
    } else {
	my $pdl_rat=$pdl_flux_fit/$pdl_model_n;
	$pdl_rat->inplace->setnantobad();
	$pdl_rat->inplace->setbadtoval(0);
	my @stats=stats($pdl_rat);
	$coeffs_0=ones(1);
	$pdl_model_0=$stats[2]*$pdl_model_n;
	$coeffs_0=$stats[2]*$coeffs_0;
#	print "@stats\n";
    }

#	my $pdl_grandom=grandom($nx);
#	my $pdl_flux_fit=$pdl_flux-$cont+$pdl_grandom*$pdl_e_flux;
#	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>1/$pdl_e_flux});
	$pdl_model_0=$pdl_model_0+$cont;
#	print "$pdl_model_0\n";
	my $n_mod_free=0;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			$a_now[$i1][1]=$coeffs_0->at($n_mod_free);
			$n_mod_free++;
		    } else {
			my $k=$link[$i1][1]-1;
			my $method=$a1[$i1][1];
			if ($method==0) {
			    $a_now[$i1][1]=$a_now[$k][1]+$a0[$i1][1];
			} else {
			    $a_now[$i1][1]=$a_now[$k][1]*$a0[$i1][1];
			}			    
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $val=$coeffs_0->at($n_mod_free);
			$a_now[$i1][$j]=$val;
			$n_mod_free++;
		    }
		}
	    }
	}
	
	
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);

	my $pdl_a_now=pdl(@a_now);


	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
	if ($chi_sq<=$chi_sq_ini) {


	    @a_out=copy_a($n_mod,\@a_now);
#	    @a_out=@a_now;
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
#	    $ii=0;
	} else {
	    if ((abs($chi_sq-$chi_sq_ini)<$d_chi_sq_goal)||($chi_sq_ini<$chi_sq_goal)) {
		$ii=$n_mc;
	    }
	}
	my $time2=get_seconds();
	my $d_time2=$time2-$time1;
#	print "$time1, $time2, $d_time2, $max_time, $ii,$n_mc\n";
	if ($d_time2>$max_time) {
	    $ii=$n_mc;
	}
	$ii++;
    }
    $chi_sq_now=$chi_sq_ini;
    my $pdl_a_out=pdl(@a_out);
#    print "$pdl_a_out\n";
#    return @a_out;
#    print "PDL MODEL ***** \n $pdl_model\n END*****\n";
    return ($chi_sq_now,$pdl_a_out,$pdl_model,$pdl_model_cont);
}




sub fit_elines_grad_rnd() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_mc=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];
    my @a_out=copy_a($n_mod,\@a);
    my $SCALE=$SCALE_INI;
    my $max_time=$_[16];

    my @dims=$pdl_flux->dims;
    my $time1=get_seconds();

    my $nx,$ny;
    if ($#dims==0) {
	$nx=$dims[0];
	$ny=1;
    }

    if ($#dims==1) {
	$nx=$dims[0];
	$ny=$dims[1];
    }

    my $i;

    my @e_stats=stats($pdl_e_flux);
    $pdl_e_flux->inplace->setvaltobad(0);
    $pdl_e_flux->inplace->setbadtoval($e_stats[0]); 



    my $pdl_model=zeroes($nx);
    my $pdl_model_cont=zeroes($nx);
    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0=zeroes($nx);
    my $pdl_model_cont_0=zeroes($nx);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
    }
    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);
    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


    my $i,$j,$ii;
    my @a_now;


    # We start the fitting loop!
    $pdl_rnd=grandom($n_mod*9*$n_mc);
    $pdl_rnd_lin=random($n_mod*9*$n_mc);
    my $ii=0;
    while ($ii<$n_mc) {

	#
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			if ($j==3) {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$rnd*($a1[$i][$j]-$a0[$i][$j])/5;
			} else {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			}
			if ($a_now[$i][$j]<$a0[$i][$j]) {
			    $a_now[$i][$j]=$a0[$i][$j];
			}
			if ($a_now[$i][$j]>$a1[$i][$j]) {
			    $a_now[$i][$j]=$a1[$i][$j];
			}
		    } else {
			$k=$link[$i][$j]-1;
			$method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
#		print "$a_now[$i][$j] $ia[$i][$j] $a0[$i][$j] $a1[$i][$j] $link[$i][$j]\n";
	    }
	}
	
	# MC test
	$pdl_model_0=zeroes($nx);
	$pdl_model_cont_0=zeroes($nx);
	for ($i=0;$i<$n_mod;$i++) {
	    my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
	    $pdl_model_0=$pdl_model_0+$pdl_tmp;
	    if ($type[$i] eq "poly1d") {
		$pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
	    }
	}
	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
	if ($chi_sq<=$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
#	    $ii=0;
	} else {
	    if ((abs($chi_sq-$chi_sq_ini)<$d_chi_sq_goal)||($chi_sq_ini<$chi_sq_goal)) {
		$ii=$n_mc;
	    }
	}
	my $time2=get_seconds();
	my $d_time2=$time2-$time1;
#	print "$time1, $time2, $d_time2, $max_time, $ii,$n_mc\n";
	if ($d_time2>$max_time) {
	    $ii=$n_mc;
	}
	$ii++;
    }
    $chi_sq_now=$chi_sq_ini;
    my $pdl_a_out=pdl(@a_out);
    return ($chi_sq_now,$pdl_a_out,$pdl_model,$pdl_model_cont);
}



sub fit_elines_TEST0() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_mc=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];
    my @a_out=copy_a($n_mod,\@a);
    my $SCALE=$SCALE_INI;
    my $max_time=$_[16];

    my @dims=$pdl_flux->dims;
    my $time1=get_seconds();

    my $nx,$ny;
    if ($#dims==0) {
	$nx=$dims[0];
	$ny=1;
    }

    if ($#dims==1) {
	$nx=$dims[0];
	$ny=$dims[1];
    }

    my $i;

    my @e_stats=stats($pdl_e_flux);
    $pdl_e_flux->inplace->setvaltobad(0);
    $pdl_e_flux->inplace->setbadtoval($e_stats[0]); 

    my $pdl_model=zeroes($nx);
    my $pdl_model_cont=zeroes($nx);
    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0,$coeffs_0;#=zeroes($nx);
    my $pdl_model_cont_0=zeroes($nx);

    my $pdl_model_n=zeroes($nx,$n_mod);

    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model_one($pdl_wave,$i,\@type,\@a);
	my $t=$pdl_model_n->slice(":,$i");
	$t .= $pdl_tmp;
    }


   ($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux,$pdl_model_n,{Weights=>1/$pdl_e_flux});
    for ($i=0;$i<$n_mod;$i++) {
	$a[$i][1]=$coeffs_0->at($i);
    }
    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);




    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


    my $i,$j,$ii;
    my @a_now;


    # We start the fitting loop!
    $pdl_rnd=grandom($n_mod*9*$n_mc);
    $pdl_rnd_lin=random($n_mod*9*$n_mc);
    my $ii=0;
    while ($ii<$n_mc) {

	#
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			if ($j==3) {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$rnd*($a1[$i][$j]-$a0[$i][$j])/5;
			} else {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			}
			if ($a_now[$i][$j]<$a0[$i][$j]) {
			    $a_now[$i][$j]=$a0[$i][$j];
			}
			if ($a_now[$i][$j]>$a1[$i][$j]) {
			    $a_now[$i][$j]=$a1[$i][$j];
			}
		    } else {
			$k=$link[$i][$j]-1;
			$method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
#		print "$a_now[$i][$j] $ia[$i][$j] $a0[$i][$j] $a1[$i][$j] $link[$i][$j]\n";
	    }
	}
	
	# MC test
#	$pdl_model_0=zeroes($nx);
#	$pdl_model_cont_0=zeroes($nx);
#	for ($i=0;$i<$n_mod;$i++) {
#	    my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
#	    $pdl_model_0=$pdl_model_0+$pdl_tmp;
#	    if ($type[$i] eq "poly1d") {
#		$pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
#	    }
#	}

	for ($i=0;$i<$n_mod;$i++) {
	    my $pdl_tmp=create_single_model_one($pdl_wave,$i,\@type,\@a_now);
	    my $t=$pdl_model_n->slice(":,$i");
	    $t .= $pdl_tmp;
	}	
	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux,$pdl_model_n,{Weights=>1/$pdl_e_flux});
	for ($i=0;$i<$n_mod;$i++) {
	    $a_now[$i][1]=$coeffs_0->at($i);
	}
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);

	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
	if ($chi_sq<=$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
#	    $ii=0;
	} else {
	    if ((abs($chi_sq-$chi_sq_ini)<$d_chi_sq_goal)||($chi_sq_ini<$chi_sq_goal)) {
		$ii=$n_mc;
	    }
	}
	my $time2=get_seconds();
	my $d_time2=$time2-$time1;
#	print "$time1, $time2, $d_time2, $max_time, $ii,$n_mc\n";
	if ($d_time2>$max_time) {
	    $ii=$n_mc;
	}
	$ii++;
    }
    $chi_sq_now=$chi_sq_ini;
    my $pdl_a_out=pdl(@a_out);
    return ($chi_sq_now,$pdl_a_out,$pdl_model,$pdl_model_cont);
}



sub fit_elines_grad_rnd_new() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_rnd=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];
    my @a_out;    
    my $SCALE=$SCALE_INI;

    my $NX;
    my @dims=$pdl_flux->dims;
    if ($#dims==0) {
	$NX=$dims[0];
#	$NY=1;
    }

    my @e_stats=stats($pdl_e_flux);
    $pdl_e_flux->inplace->setvaltobad(0);
    $pdl_e_flux->inplace->setbadtoval($e_stats[0]); 
 

    my $i;

    my $pdl_model=zeroes($NX);
    my $pdl_model_cont=zeroes($NX);
    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0=zeroes($NX);
    my $pdl_model_cont_0=zeroes($NX);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
    }
    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);
    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=3*($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


    my $i,$j,$ii;
    my @a_now;


    $pdl_rnd=grandom($n_mod*9*$n_mc);
    $pdl_rnd_lin=0.8+0.4*random($n_mod*9*$n_mc);





#
# 1st we derive the redshift!
#
    $new_mc=int($n_mc/2);
    for ($ii=0;$ii<$new_mc;$ii++) {
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			if ($type[$i] ne "poly1d") {
			    if ($j==3) {
				$a_now[$i][$j]=$a0[$i][$j]+$rnd_lin*($a1[$i][$j]-$a0[$i][$j])*$ii/$new_mc;
			    } else {
				$a_now[$i][$j]=$a_out[$i][$j];#+$SCALE*$a_out[$i][$j]*$rnd;
			    }
			} else {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			}

			if ($a_now[$i][$j]<$a0[$i][$j]) {
			    $a_now[$i][$j]=$a0[$i][$j];
			}
			if ($a_now[$i][$j]>$a1[$i][$j]) {
			    $a_now[$i][$j]=$a1[$i][$j];
			}

		    } else {
			my $k=$link[$i][$j]-1;
			my $method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
	    }
	}

	my $n_mod_free=0;
	my $i1,$j1;
	my $pdl_model_n;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			if ($n_mod_free==0) {
			    $pdl_model_n=$pdl_tmp;
			} else {
			    $pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			}
			$n_mod_free++;
		    } else {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			my $ii=$link[$i1][1]-1;
			my $t=$pdl_model_n->slice(":,$ii");
			$pdl_tmp=$pdl_tmp*$a0[$i1][1];		
			$t .= $t +$pdl_tmp;		   
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {	    
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
			$pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			$n_mod_free++;
		    }
		}
	    }
	}


	my $pdl_grandom=grandom($NX);
	my $pdl_flux_fit=$pdl_flux-$cont+$pdl_grandom*$pdl_e_flux;
#	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>1/$pdl_e_flux});

    @dim_n=$pdl_model_n->dims();
    if ($dim_n[1]>1) {
#	print "$pdl_flux_fit\n $pdl_model_n\n";

#	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>1/$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n);#,{Weights=>1/$pdl_e_flux});
#	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n);#,{Weights=>1/$pdl_e_flux});
    } else {
	my $pdl_rat=$pdl_flux_fit/$pdl_model_n;
	$pdl_rat->inplace->setnantobad();
	$pdl_rat->inplace->setbadtoval(0);
	my @stats=stats($pdl_rat);
	$coeffs_0=ones(1);
	$pdl_model_0=$stats[2]*$pdl_model_n;
	$coeffs_0=$stats[2]*$coeffs_0;
    }


	$pdl_model_0=$pdl_model_0+$cont;
#	print "$pdl_model_0\n";
	my $n_mod_free=0;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			$a_now[$i1][1]=$coeffs_0->at($n_mod_free);

			if ($a_now[$i1][1]<$a0[$i1][1]) {
			    $a_now[$i1][1]=$a0[$i1][1];
			}
			if ($a_now[$i1][1]>$a1[$i1][1]) {
			    $a_now[$i1][1]=$a1[$i1][1];
			}
			$n_mod_free++;
		    } else {
			my $k=$link[$i1][1]-1;
			my $method=$a1[$i1][1];
			if ($method==0) {
			    $a_now[$i1][1]=$a_now[$k][1]+$a0[$i1][1];
			} else {
			    $a_now[$i1][1]=$a_now[$k][1]*$a0[$i1][1];
			}			    
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $val=$coeffs_0->at($n_mod_free);
			$a_now[$i1][$j]=$val;
			$n_mod_free++;
		    }
		}
	    }
	}


	my $pdl_model_0=zeroes($NX);
	my $pdl_model_cont_0=zeroes($NX);
	for ($i=0;$i<$n_mod;$i++) {
	    my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
	    $pdl_model_0=$pdl_model_0+$pdl_tmp;
	}
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);

#	plot_results(1,$pdl_wave,pdl($pdl_flux,$pdl_model,$pdl_model_0,$pdl_e_flux),"junk","redshift");	
	

	my $pdl_a_now=pdl(@a_now);

	
	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
#	print "$ii/$new_mc $a_now[0][3] ($chi_sq<=$chi_sq_ini)\n";
	if ($chi_sq<=$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
#	    @a_out=@a_now;
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
	}
#	my $time2=get_seconds();
#	my $d_time2=$time2-$time1;
#	print "$time1, $time2, $d_time2, $max_time, $ii,$n_mc\n";
#	if ($d_time2>$max_time) {
#	    $ii=$n_mc;
#	}

#    plot_results(1,$pdl_wave,pdl($pdl_flux,$pdl_model,$pdl_model_0,$pdl_e_flux),"junk","redshift");

#	print "$ii/$new_mc $chi_sq<$chi_sq_ini $n_points,$n_mod\n";#$pdl_e_flux\n";

#	for ($h=0;$h<$NX;$h++) {
#	    my $val=$pdl_e_flux->at($h);
#	    if ($val==0) {
#		print "$val\n";
#	    }
#	}


   }

#    print_a_val($n_mod,\@a0,\@a1,\@test);

    $ii=0;

    $chi_sq_ini=1e12;

    @a_now=copy_a($n_mod,\@a_out);


#    print "-------------\n";
#   print_a_val($n_mod,\@a_out,\@a_now,\@type,$chi_sq_now);

#
# 2nd we derive the sigma
#
    $new_mc=int($n_mc/3);
    for ($ii=0;$ii<$new_mc;$ii++) {
	#
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			if ($type[$i] ne "poly1d") {
			    if ($j==2) {
				$a_now[$i][$j]=$a0[$i][$j]+$rnd_lin*($a1[$i][$j]-$a0[$i][$j])*$ii/$new_mc;#$a_out[$i][$j];#+$SCALE*$a_out[$i][$j]*$rnd;
			    } else {
				$a_now[$i][$j]=$a_out[$i][$j];
			    }
			} else {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			}

			if ($a_now[$i][$j]<$a0[$i][$j]) {
			    $a_now[$i][$j]=$a0[$i][$j];
			}
			if ($a_now[$i][$j]>$a1[$i][$j]) {
			    $a_now[$i][$j]=$a1[$i][$j];
			}

		    } else {
#		    print "PASO $i,$j $link[$i][$j]\n";

			my $k=$link[$i][$j]-1;
			my $method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
	    }
	}

	my $n_mod_free=0;
	my $i1,$j1;
	my $pdl_model_n;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			if ($n_mod_free==0) {
			    $pdl_model_n=$pdl_tmp;
			} else {
			    $pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			}
			$n_mod_free++;
		    } else {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			my $ii=$link[$i1][1]-1;
			my $t=$pdl_model_n->slice(":,$ii");
			$pdl_tmp=$pdl_tmp*$a0[$i1][1];		
			$t .= $t +$pdl_tmp;		   
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {	    
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
			$pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			$n_mod_free++;
		    }
		}
	    }
	}


	my $pdl_grandom=grandom($NX);
	my $pdl_flux_fit=$pdl_flux-$cont+$pdl_grandom*$pdl_e_flux;

    @dim_n=$pdl_model_n->dims();
    if ($dim_n[1]>1) {
	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n);#,{Weights=>1/$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
    } else {
	my $pdl_rat=$pdl_flux_fit/$pdl_model_n;
	$pdl_rat->inplace->setnantobad();
	$pdl_rat->inplace->setbadtoval(0);
	my @stats=stats($pdl_rat);
	$coeffs_0=ones(1);
	$pdl_model_0=$stats[2]*$pdl_model_n;
	$coeffs_0=$stats[2]*$coeffs_0;
    }


	$pdl_model_0=$pdl_model_0+$cont;
#	print "$pdl_model_0\n";
	my $n_mod_free=0;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
#		    print "** PASO $i,$j $link[$i][$j]\n";
		    if ($link[$i1][1]==-1) {
			$a_now[$i1][1]=$coeffs_0->at($n_mod_free);
			$n_mod_free++;
		    } else {
			my $k=$link[$i1][1]-1;
			my $method=$a1[$i1][1];
			if ($method==0) {
			    $a_now[$i1][1]=$a_now[$k][1]+$a0[$i1][1];
			} else {
			    $a_now[$i1][1]=$a_now[$k][1]*$a0[$i1][1];
			}			    
#			print "** $a_now[$i1][1]=$a_now[$k][1]*$a0[$i1][1] \n"; 
		    }

		}

#
# Check they are in the range!
#
		for ($j=0;$j<9;$j++) {
		    if (($ia[$i1][$j]==1)&&($link[$i1][$j]==-1)) {
			if ($a_now[$i1][$j]<$a0[$i1][$j]) {
			    $a_now[$i1][$j]=$a0[$i1][$j];
			}
			if ($a_now[$i1][$j]>$a1[$i1][$j]) {
			    $a_now[$i1][$j]=$a1[$i1][$j];
			}
		    }
		}

	    }
	    if ($type[$i1] eq "poly1d") {
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $val=$coeffs_0->at($n_mod_free);
			$a_now[$i1][$j]=$val;
			$n_mod_free++;
		    }
		}
	    }
	}
	
	
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);

	my $pdl_a_now=pdl(@a_now);


	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
#	print "$ii/$new_mc $a_now[0][4] ($chi_sq<=$chi_sq_ini)\n";
	if ($chi_sq<=$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
#	    @a_out=@a_now;
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
	}
    }


    $ii=0;

#    $chi_sq_ini=1e12;

    @a_now=copy_a($n_mod,\@a_out);
    my @a_now_lin=copy_a($n_mod,\@a_now);

    # We start the fitting loop!

    my $ii=0;
    my $SCALE_IN=$SCALE;
    $new_mc=0;
#    $ii=$n_mc;
    while ($ii<$new_mc) {

	#
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			my $A1=$a1[$i][$j];
			my $A0=$a0[$i][$j];
			if ($a1[$i][$j]>1.3*$a_out[$i][$j]) {
			    $A1=1.3*$a_out[$i][$j];
			}
			if ($a0[$i][$j]<0.7*$a_out[$i][$j]) {
			    $A0=0.7*$a_out[$i][$j];
			}

#			print "$A1,$A0\n";

			if ($type[$i] eq "eline") {
			    if ($j==3) {
				$a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$rnd*($A1-$A0)/(5*$new_mc);
			    } else {
				$a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			    }

			if ($a_now[$i][$j]<$A0) {
			    $a_now[$i][$j]=$A0;
			}
			if ($a_now[$i][$j]>$A1) {
			    $a_now[$i][$j]=$A1;
			}



			} else {
# No variation!
			    $a_now[$i][$j]=$a_out[$i][$j];#+$SCALE*$a_out[$i][$j]*$rnd*0.0001;
			}
			
			

		    } else {
			$k=$link[$i][$j]-1;
			$method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
	    }
#	    print "$type[$i] $a_now[$i][0]\n";
	}

	# MC test
	$pdl_model_0=zeroes($NX);
	$pdl_model_cont_0=zeroes($NX);
	for ($i=0;$i<$n_mod;$i++) {
	    my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
	    $pdl_model_0=$pdl_model_0+$pdl_tmp;
	}
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);
	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
	if ($chi_sq<$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
	    $SCALE=$SCALE*0.99;
	    if ($SCALE<0.1*$SCALE_IN) {
		$SCALE=$SCALE_IN*0.1;
	    }
	} else {
	    $SCALE=$SCALE_IN;

	    if ((abs($chi_sq-$chi_sq_ini)<$d_chi_sq_goal)||($chi_sq_ini<$chi_sq_goal)) {
		$ii=$n_mc;
	    }
	}
	$ii++;
    }
    $chi_sq_now=$chi_sq_ini;

#
# We force lineal!
#
    my @a_out=copy_a($n_mod,\@a_now_lin);
#    print "@a_out\n";
#    print "*****\n";
#    print_a_val($n_mod,\@a_out,\@a_now_lin,\@type,$chi_sq_now);
    my $pdl_a_out=pdl(@a_out);
#    print "*****\n";
    return ($chi_sq_now,$pdl_a_out,$pdl_model,$pdl_model_cont);
#    return @a_out;
}

sub fit_elines_grad_rnd_new_guided() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_rnd=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];
    my @a_out;    
    my $SCALE=$SCALE_INI;
    my $g_v=$_[16];
    my $g_d=$_[17];

    my $NX;
    my @dims=$pdl_flux->dims;
    if ($#dims==0) {
	$NX=$dims[0];
#	$NY=1;
    }

    my @e_stats=stats($pdl_e_flux);
    $pdl_e_flux->inplace->setvaltobad(0);
    $pdl_e_flux->inplace->setbadtoval($e_stats[0]); 
 

    my $i;

    my $pdl_model=zeroes($NX);
    my $pdl_model_cont=zeroes($NX);
    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0=zeroes($NX);
    my $pdl_model_cont_0=zeroes($NX);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
    }
    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);
    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


    my $i,$j,$ii;
    my @a_now;


    $pdl_rnd=grandom($n_mod*9*$n_mc);
    $pdl_rnd_lin=0.8+0.4*random($n_mod*9*$n_mc);





#
# 1st we derive the redshift!
#
    if ($g_v==0) {
	$new_mc=3;
    } else {
	$new_mc=int($n_mc/2);
    }
    for ($ii=0;$ii<$new_mc;$ii++) {
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			if ($type[$i] ne "poly1d") {
			    if ($j==3) {
				$a_now[$i][$j]=$a0[$i][$j]+$rnd_lin*($a1[$i][$j]-$a0[$i][$j])*$ii/$new_mc;
			    } else {
				$a_now[$i][$j]=$a_out[$i][$j];#+$SCALE*$a_out[$i][$j]*$rnd;
			    }
			} else {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			}

			if ($a_now[$i][$j]<$a0[$i][$j]) {
			    $a_now[$i][$j]=$a0[$i][$j];
			}
			if ($a_now[$i][$j]>$a1[$i][$j]) {
			    $a_now[$i][$j]=$a1[$i][$j];
			}

		    } else {
			my $k=$link[$i][$j]-1;
			my $method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
	    }
	}

	my $n_mod_free=0;
	my $i1,$j1;
	my $pdl_model_n;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			if ($n_mod_free==0) {
			    $pdl_model_n=$pdl_tmp;
			} else {
			    $pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			}
			$n_mod_free++;
		    } else {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			my $ii=$link[$i1][1]-1;
			my $t=$pdl_model_n->slice(":,$ii");
			$pdl_tmp=$pdl_tmp*$a0[$i1][1];		
			$t .= $t +$pdl_tmp;		   
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {	    
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
			$pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			$n_mod_free++;
		    }
		}
	    }
	}


	my $pdl_grandom=grandom($NX);
	my $pdl_flux_fit=$pdl_flux-$cont;#+$pdl_grandom*$pdl_e_flux;
#	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>1/$pdl_e_flux});

    @dim_n=$pdl_model_n->dims();
    if ($dim_n[1]>1) {
#	print "$pdl_flux_fit\n $pdl_model_n\n";
	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n,{Weights=>1/$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
    } else {
	my $pdl_rat=$pdl_flux_fit/$pdl_model_n;
	$pdl_rat->inplace->setnantobad();
	$pdl_rat->inplace->setbadtoval(0);
	my @stats=stats($pdl_rat);
	$coeffs_0=ones(1);
	$pdl_model_0=$stats[2]*$pdl_model_n;
	$coeffs_0=$stats[2]*$coeffs_0;
    }


	$pdl_model_0=$pdl_model_0+$cont;
#	print "$pdl_model_0\n";
	my $n_mod_free=0;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			$a_now[$i1][1]=$coeffs_0->at($n_mod_free);

			if ($a_now[$i1][1]<$a0[$i1][1]) {
			    $a_now[$i1][1]=$a0[$i1][1];
			}
			if ($a_now[$i1][1]>$a1[$i1][1]) {
			    $a_now[$i1][1]=$a1[$i1][1];
			}


			$n_mod_free++;
		    } else {
			my $k=$link[$i1][1]-1;
			my $method=$a1[$i1][1];
			if ($method==0) {
			    $a_now[$i1][1]=$a_now[$k][1]+$a0[$i1][1];
			} else {
			    $a_now[$i1][1]=$a_now[$k][1]*$a0[$i1][1];
			}			    
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $val=$coeffs_0->at($n_mod_free);
			$a_now[$i1][$j]=$val;
			$n_mod_free++;
		    }
		}
	    }
	}


	my $pdl_model_0=zeroes($NX);
	my $pdl_model_cont_0=zeroes($NX);
	for ($i=0;$i<$n_mod;$i++) {
	    my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
	    $pdl_model_0=$pdl_model_0+$pdl_tmp;
	}
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);

#	plot_results(1,$pdl_wave,pdl($pdl_flux,$pdl_model,$pdl_model_0,$pdl_e_flux),"junk","redshift");	
	

	my $pdl_a_now=pdl(@a_now);

	
	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
#	print "$ii/$new_mc $a_now[0][3] ($chi_sq<=$chi_sq_ini)\n";
	if ($chi_sq<=$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
#	    @a_out=@a_now;
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
	}
#	my $time2=get_seconds();
#	my $d_time2=$time2-$time1;
#	print "$time1, $time2, $d_time2, $max_time, $ii,$n_mc\n";
#	if ($d_time2>$max_time) {
#	    $ii=$n_mc;
#	}

#    plot_results(1,$pdl_wave,pdl($pdl_flux,$pdl_model,$pdl_model_0,$pdl_e_flux),"junk","redshift");

#	print "$ii/$new_mc $chi_sq<$chi_sq_ini $n_points,$n_mod\n";#$pdl_e_flux\n";

#	for ($h=0;$h<$NX;$h++) {
#	    my $val=$pdl_e_flux->at($h);
#	    if ($val==0) {
#		print "$val\n";
#	    }
#	}


   }

#    print_a_val($n_mod,\@a0,\@a1,\@test);

    $ii=0;

    $chi_sq_ini=1e12;

    @a_now=copy_a($n_mod,\@a_out);

#
# 2nd we derive the sigma
#
    if ($g_d==0) {
	$new_mc=3;
    } else {
	$new_mc=int($n_mc/3);
    }
#    $new_mc=int($n_mc/3);
    for ($ii=0;$ii<$new_mc;$ii++) {
	#
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			if ($type[$i] ne "poly1d") {
			    if ($j==2) {
				$a_now[$i][$j]=$a0[$i][$j]+$rnd_lin*($a1[$i][$j]-$a0[$i][$j])*$ii/$new_mc;#$a_out[$i][$j];#+$SCALE*$a_out[$i][$j]*$rnd;
			    } else {
				$a_now[$i][$j]=$a_out[$i][$j];
			    }
			} else {
			    $a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			}

			if ($a_now[$i][$j]<$a0[$i][$j]) {
			    $a_now[$i][$j]=$a0[$i][$j];
			}
			if ($a_now[$i][$j]>$a1[$i][$j]) {
			    $a_now[$i][$j]=$a1[$i][$j];
			}

		    } else {
#			print "PASO $i,$j\n";
			my $k=$link[$i][$j]-1;
			my $method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    

		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
	    }
	}

	my $n_mod_free=0;
	my $i1,$j1;
	my $pdl_model_n;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			if ($n_mod_free==0) {
			    $pdl_model_n=$pdl_tmp;
			} else {
			    $pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			}
			$n_mod_free++;
		    } else {
			my $pdl_tmp=create_single_model_one($pdl_wave,$i1,\@type,\@a_now);
			my $ii=$link[$i1][1]-1;
			my $t=$pdl_model_n->slice(":,$ii");
			$pdl_tmp=$pdl_tmp*$a0[$i1][1];		
			$t .= $t +$pdl_tmp;		   
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {	    
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $pdl_tmp=create_single_model_poly($pdl_wave,$j);
			$pdl_model_n=$pdl_model_n->glue(1,$pdl_tmp);
			$n_mod_free++;
		    }
		}
	    }
	}


	my $pdl_grandom=grandom($NX);
	my $pdl_flux_fit=$pdl_flux-$cont;#+$pdl_grandom*$pdl_e_flux;

    @dim_n=$pdl_model_n->dims();
    if ($dim_n[1]>1) {
	($pdl_model_0, $coeffs_0) = linfit1d($pdl_flux_fit,$pdl_model_n);#,{Weights=>1/$pdl_e_flux});#,{Weights=>1/$pdl_e_flux});
    } else {
	my $pdl_rat=$pdl_flux_fit/$pdl_model_n;
	$pdl_rat->inplace->setnantobad();
	$pdl_rat->inplace->setbadtoval(0);
	my @stats=stats($pdl_rat);
	$coeffs_0=ones(1);
	$pdl_model_0=$stats[2]*$pdl_model_n;
	$coeffs_0=$stats[2]*$coeffs_0;
    }


	$pdl_model_0=$pdl_model_0+$cont;
#	print "$pdl_model_0\n";
	my $n_mod_free=0;
	for ($i1=0;$i1<$n_mod;$i1++) {
	    if ($type[$i1] eq "eline") {
		if ($ia[$i1][1]==1) {
		    if ($link[$i1][1]==-1) {
			$a_now[$i1][1]=$coeffs_0->at($n_mod_free);
			$n_mod_free++;
		    } else {
			my $k=$link[$i1][1]-1;
			my $method=$a1[$i1][1];
			if ($method==0) {
			    $a_now[$i1][1]=$a_now[$k][1]+$a0[$i1][1];
			} else {
			    $a_now[$i1][1]=$a_now[$k][1]*$a0[$i1][1];
			}			    
		    }
		}
	    }
	    if ($type[$i1] eq "poly1d") {
		for ($j=0;$j<9;$j++) {
		    if ($ia[$i1][$j]==1) {
			my $val=$coeffs_0->at($n_mod_free);
			$a_now[$i1][$j]=$val;
			$n_mod_free++;
		    }
		}
	    }
	}
	
	
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);

	my $pdl_a_now=pdl(@a_now);


	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
#	print "$ii/$new_mc $a_now[0][4] ($chi_sq<=$chi_sq_ini)\n";
	if ($chi_sq<=$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
#	    @a_out=@a_now;
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
	}
    }


    $ii=0;

    $chi_sq_ini=1e12;

    @a_now=copy_a($n_mod,\@a_out);



    # We start the fitting loop!

    $new_mc=int($n_mc/3);
 
    my $ii=0;
    my $SCALE_IN=$SCALE;
#    $ii=$n_mc;
    while ($ii<$new_mc) {

	#
	# We change slighly the parameters
	#
	for ($i=0;$i<$n_mod;$i++) {
	    for ($j=0;$j<9;$j++) {
		$rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
		$rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
		if ($ia[$i][$j]==1) {
		    if ($link[$i][$j]==-1) {
			my $A1=$a1[$i][$j];
			my $A0=$a0[$i][$j];
			if ($a1[$i][$j]>1.3*$a_out[$i][$j]) {
			    $A1=1.3*$a_out[$i][$j];
			}
			if ($a0[$i][$j]<0.7*$a_out[$i][$j]) {
			    $A0=0.7*$a_out[$i][$j];
			}

#			print "$A1,$A0\n";

			if ($type[$i] eq "eline") {
			    if ($j==3) {
				$a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$rnd*($A1-$A0)/(5*$new_mc);
			    } else {
				$a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
			    }

			if ($a_now[$i][$j]<$A0) {
			    $a_now[$i][$j]=$A0;
			}
			if ($a_now[$i][$j]>$A1) {
			    $a_now[$i][$j]=$A1;
			}



			} else {
# No variation!
			    $a_now[$i][$j]=$a_out[$i][$j];#+$SCALE*$a_out[$i][$j]*$rnd*0.0001;
			}
			
			

		    } else {
			$k=$link[$i][$j]-1;
			$method=$a1[$i][$j];
			if ($method==0) {
			    $a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
			} else {
			    $a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
			}			    
		    }
		} else {
		    $a_now[$i][$j]=$a_out[$i][$j];
		}	    
	    }
#	    print "$type[$i] $a_now[$i][0]\n";
	}

	# MC test
	$pdl_model_0=zeroes($NX);
	$pdl_model_cont_0=zeroes($NX);
	for ($i=0;$i<$n_mod;$i++) {
	    my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
	    $pdl_model_0=$pdl_model_0+$pdl_tmp;
	}
	$pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);
	my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
	my $chi_sq=sum($pdl_chi_now);
	$chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
	if ($chi_sq<$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
	    $SCALE=$SCALE*0.99;
	    if ($SCALE<0.1*$SCALE_IN) {
		$SCALE=$SCALE_IN*0.1;
	    }
	} else {
	    $SCALE=$SCALE_IN;

	    if ((abs($chi_sq-$chi_sq_ini)<$d_chi_sq_goal)||($chi_sq_ini<$chi_sq_goal)) {
		$ii=$n_mc;
	    }
	}
	$ii++;
    }
    $chi_sq_now=$chi_sq_ini;

    my $pdl_a_out=pdl(@a_out);
    return ($chi_sq_now,$pdl_a_out,$pdl_model,$pdl_model_cont);
#    return @a_out;
}






sub fit_elines_grad_rnd_NEW() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_mc=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];    
    my $SCALE=$SCALE_INI;

    my $nx,$ny;
    my @dims=$pdl_flux->dims;
    if ($#dims==0) {
	$nx=$dims[0];
	$ny=1;
    }

    if ($#dims==1) {
	$nx=$dims[0];
	$ny=$dims[1];
    }

    my $i;


    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0=zeroes($nx);
    my $pdl_model_cont_0=zeroes($nx);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
	if ($type[$i] eq "poly1d") {
	    $pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
	}
    }
#    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);
    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


    my $i,$j,$ii;
    my @a_now;

    my $pdl_rnd=grandom($n_mod*9*$n_mc);
    my $pdl_rnd_lin=random($n_mod*9*$n_mc);

    # We start the fitting loop!
    my $ii=0;
    while ($ii<$n_mc) {
	my $pdl_chi_now;
	my $chi_sq;
	my $II;
	my @childs;
	my $pid;


	my $junk=fit_elines_rnd_single($pdl_wave,$pdl_flux,$pdl_e_flux,$n_mod,$chi_goal,$d_chi_goal,\@type,\@a_out,\@ia,\@a0,\@a1,\@link,$n_mc,$pdl_masked,$def,$scale_ini,$pdl_rnd,$pdl_rnd_lin,$ii);
	my @a_now=@$junk;
	$chi_sq=$a_now[$n_mod][0];

#	$chi_sq=$chi_single;

#	print "$ii,$n_mc,$chi_sq,$chi_sq_ini\n";

	
	if ($chi_sq<$chi_sq_ini) {
	    @a_out=copy_a($n_mod,\@a_now);
	    $chi_sq_ini=$chi_sq;
	    $pdl_model=$pdl_model_0;
	    $pdl_model_cont=$pdl_model_cont_0;
	} else {
	    if ((abs($chi_sq-$chi_sq_ini)<$d_chi_sq_goal)||($chi_sq_ini<$chi_sq_goal)) {
		$ii=$n_mc;
	    }
	}

	$ii++;
    }
    $chi_sq_now=$chi_sq_ini;
    return @a_out;
}



sub fit_elines_grad_rnd_fork() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_mc=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];
    my $NC=$_[16];
    my @a_out;    
    my $SCALE=$SCALE_INI;

    my @a_fork;
    
    my $i; 

    my $nx,$ny;
    my @dims=$pdl_flux->dims;
    if ($#dims==0) {
	$nx=$dims[0];
	$ny=1;
    }

    if ($#dims==1) {
	$nx=$dims[0];
	$ny=$dims[1];
    }

    my $i;


    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0=zeroes($nx);
    my $pdl_model_cont_0=zeroes($nx);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
	if ($type[$i] eq "poly1d") {
	    $pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
	}
    }
#    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);
    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


    my $i,$j,$ii;
    my @a_now;

    my $pdl_rnd=grandom($n_mod*9*$n_mc);
    my $pdl_rnd_lin=random($n_mod*9*$n_mc);

    # We start the fitting loop!
    my $ii=0;
    my @a_now_final;
    my @a_now_single;

    while ($ii<$n_mc) {
	my $pdl_chi_now;
	my $chi_sq;
	my $II;
	my @childs;
	my $pid;

# $forks[$i] = async { return sum($i,$sum); };

	for ($II=0;$II<$NC;$II++) {
	    $a_now_single[$II]= async {     my $pdl_rnd=grandom($n_mod*9*$n_mc);
					    my $pdl_rnd_lin=random($n_mod*9*$n_mc);
					    fit_elines_rnd_single($pdl_wave,$pdl_flux,$pdl_e_flux,$n_mod,$chi_goal,$d_chi_goal,\@type,\@a_out,\@ia,\@a0,\@a1,\@link,$n_mc,$pdl_masked,$def,$scale_ini,$pdl_rnd,$pdl_rnd_lin,$ii);
	    }
	}
	

	for ($II=0;$II<$NC;$II++) {
	    my $junk = $a_now_single[$II]->join();
	    @a_now=@$junk;
	    $chi_sq=$a_now[$n_mod][0];
#	    print "$ii,$II $chi_sq\n";
	    if ($chi_sq<$chi_sq_ini) {
		@a_out=copy_a($n_mod,\@a_now);
		$chi_sq_ini=$chi_sq;
		$pdl_model=$pdl_model_0;
		$pdl_model_cont=$pdl_model_cont_0;
	    } else {
		if ((abs($chi_sq-$chi_sq_ini)<$d_chi_sq_goal)||($chi_sq_ini<$chi_sq_goal)) {
		    $ii=$n_mc;
		}
	    }
	    
	}
	$ii=$ii+$NC;
	$chi_sq_now=$chi_sq_ini;
    }
    return @a_out;
}



sub fit_elines_rnd_single() {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_e_flux=$_[2];
    my $n_mod=$_[3];
    my $chi_goal=$_[4];
    my $d_chi_goal=$_[5];
    my $junk=$_[6];
    my @type=@$junk;
    my $junk=$_[7];
    my @a=@$junk;
    my $junk=$_[8];
    my @ia=@$junk;
    my $junk=$_[9];
    my @a0=@$junk;
    my $junk=$_[10];
    my @a1=@$junk;
    my $junk=$_[11];
    my @link=@$junk;
    my $n_mc=$_[12];
    my $pdl_masked=$_[13];
    my $def=$_[14];
    my $SCALE_INI=$_[15];
    my @a_out;    
    my $SCALE=$SCALE_INI;
    my $pdl_rnd=$_[16];
    my $pdl_rnd_lin=$_[17];
    my $ii=$_[18];
    my @dims=$pdl_flux->dims;
    my $nx,$ny;
    if ($#dims==0) {
	$nx=$dims[0];
	$ny=1;
    }

    if ($#dims==1) {
	$nx=$dims[0];
	$ny=$dims[1];
    }

    my $i;


    # First Guess
    my $n_points=sum($pdl_masked);
    my $pdl_model_0=zeroes($nx);
    my $pdl_model_cont_0=zeroes($nx);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
	if ($type[$i] eq "poly1d") {
	    $pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
	}
    }
#    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a);
    my $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    my $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_sq_ini=$chi_sq;
    @a_out=copy_a($n_mod,\@a);


    my $i,$j,$ii;
    my @a_now;


    # We start the fitting loop!

    my $pdl_chi_now;
    my $chi_sq;
    my $II;
    my @childs;
    my $pid;


#    print_a_final($n_mod,\@a_out,\@type,$chi_sq_now);

    
	# We change slighly the parameters
    #
    for ($i=0;$i<$n_mod;$i++) {
	for ($j=0;$j<9;$j++) {
	    $rnd=$pdl_rnd->at(($j+(9*$i))*$ii);
	    $rnd_lin=$pdl_rnd_lin->at(($j+(9*$i))*$ii);
	    if ($ia[$i][$j]==1) {
		if ($link[$i][$j]==-1) {
		    if ($j==3) {
			$a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$rnd*($a1[$i][$j]-$a0[$i][$j])/5;
		    } else {
			$a_now[$i][$j]=$a_out[$i][$j]+$SCALE*$a_out[$i][$j]*$rnd;
		    }
		    if ($a_now[$i][$j]<$a0[$i][$j]) {
			$a_now[$i][$j]=$a0[$i][$j];
		    }
		    if ($a_now[$i][$j]>$a1[$i][$j]) {
			$a_now[$i][$j]=$a1[$i][$j];
		    }
		} else {
		    $k=$link[$i][$j]-1;
		    $method=$a1[$i][$j];
		    if ($method==0) {
			$a_now[$i][$j]=$a_now[$k][$j]+$a0[$i][$j];
		    } else {
			$a_now[$i][$j]=$a_now[$k][$j]*$a0[$i][$j];
		    }			    
		}
	    } else {
		$a_now[$i][$j]=$a_out[$i][$j];
	    }	    
#	    print "$i,$j $a_now[$i][$j]\n";
	}
    }
    
    # MC test
    $pdl_model_0=zeroes($nx);
    $pdl_model_cont_0=zeroes($nx);
    for ($i=0;$i<$n_mod;$i++) {
	my $pdl_tmp=create_single_model($pdl_wave,$i,\@type,\@a_now);
	$pdl_model_0=$pdl_model_0+$pdl_tmp;
	if ($type[$i] eq "poly1d") {
	    $pdl_model_cont_0=$pdl_model_cont_0+$pdl_tmp;
	}
    }
#    $pdl_model_cont_0=create_single_model($pdl_wave,($n_mod-1),\@type,\@a_now);
    $pdl_chi_now=$pdl_masked*(($pdl_flux-$pdl_model_0)**2)/($pdl_e_flux**2);
    $chi_sq=sum($pdl_chi_now);
    $chi_sq=($chi_sq/($n_points-$n_mod-1))**0.5;
    $chi_single=$chi_sq;
    $a_now[$n_mod][0]=$chi_sq;
#    print "@a_now\n";
    #$chi_sq_now=$chi_sq;
    #print "$chi_sq_now\n";

    return \@a_now;
}

sub plot_results_old {
    my $plot=$_[0];
    my $pdl_wave=$_[1];
    my $pdl_output=$_[2];
    my $output_name=$_[3];
    my $title=$_[4];
    my $min_wave,$max_wave;
    my $dev_plot;
    my $y_min,$y_max;
    my @stats;
    my $nx,$ny;
    my @wave_now=list($pdl_wave);
    my $i,$j;
    my $nx_i,$ny_i;

    if ($plot>0) {
	$dev_plot="/null";
	if ($plot==1) {
	    $dev_plot="/xs";
	} else {
	    if ($plot==3) {
		$dev_plot=$output_name.".png/PNG";
	    } else {
		$dev_plot=$output_name.".ps/CPS";
	    }
	}
#	print "PLOT=$plot,DEV=$dev_plot\n";
	@stats=stats($pdl_output->slice(":,(0)"));
	$y_min=-0.5*$stats[0]-0.75*$stats[1];
	$y_max=2*$stats[0]+6*$stats[1];
	($nx_i,$ny_i)=$pdl_output->dims();
	$min_wave=$wave_now[0];
	$max_wave=$wave_now[$nx-1];

	pgbegin(0,$dev_plot,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($min_wave,$max_wave,$y_min,$y_max,0,0);
	pgsch(0.9);           # Set character height
	pglabel("Wavelength","Flux","$title");
	pgsch(0.5);
	for ($j=0;$j<$ny_i;$j++) {
	    my $t=$pdl_output->slice(":,($j)");
	    pgsci(1+$j);
	    my @flux=list($t);
#	    print "$j/$ny @flux @wave_now\n";
	    pgline($nx,\@wave_now,\@flux);
#	    print "plot $j\n"; <stdin>;
	}
	pgsci(1);
	pgclose;
	pgend;
    }
}


sub plot_results_min_max {
    my $plot=$_[0];
    my $pdl_wave=$_[1];
    my $pdl_output=$_[2];
    my $output_name=$_[3];
    my $title=$_[4];
    my $y_min=$_[5];
    my $y_max=$_[6];

    my $min_wave,$max_wave;
    my $dev_plot;
    my @stats;
    my $nx,$ny;
    my @wave_now=list($pdl_wave);
    my $i,$j;

    if ($plot>0) {
	$dev_plot="/null";
	if ($plot==1) {
	    $dev_plot="/xs";
	} else {
	    if ($plot==3) {
		$dev_plot=$output_name.".png/PNG";
	    } else {
		$dev_plot=$output_name.".ps/CPS";
	    }
	}
#	print "PLOT=$plot,DEV=$dev_plot\n";
	@stats=stats($pdl_output->slice(":,(0)"));
	#$y_min=-0.5*$stats[0]-0.75*$stats[1];
	#$y_max=2*$stats[0]+6*$stats[1];
	($nx,$ny)=$pdl_output->dims();
	$min_wave=$wave_now[0];
	$max_wave=$wave_now[$nx-1];

	pgbegin(0,$dev_plot,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($min_wave,$max_wave,$y_min,$y_max,0,0);
	pgsch(0.9);           # Set character height
	pglabel("Wavelength","Flux","$title");
	pgsch(0.5);
	for ($j=0;$j<$ny;$j++) {
	    my $t=$pdl_output->slice(":,($j)");
	    pgsci(1+$j);
	    my @flux=list($t);
	    pgline($nx,\@wave_now,\@flux);
	}
	pgsci(1);
	pgclose;
	pgend;
    }
}


sub plot_results {
    my $plot=$_[0];
    my $pdl_wave=$_[1];
    my $pdl_output=$_[2];
    my $output_name=$_[3];
    my $title=$_[4];
    my $y_min=$_[5];
    my $y_max=$_[6];

    my $min_wave,$max_wave;
    my $dev_plot;
    my @stats;
    my $nx,$ny;
    my @wave_now=list($pdl_wave);
    my $i,$j;

    if ($plot>0) {
	$dev_plot="/null";
	if ($plot==1) {
	    $dev_plot="/xs";
	} else {
	    if ($plot==3) {
		$dev_plot=$output_name.".png/PNG";
	    } else {
		$dev_plot=$output_name.".ps/CPS";
	    }
	}
#	print "PLOT=$plot,DEV=$dev_plot\n";
	@stats=stats($pdl_output->slice(":,(0)"));
	$y_min=-0.5*$stats[0]-0.75*$stats[1];
	$y_max=2*$stats[0]+6*$stats[1];
	($nx,$ny)=$pdl_output->dims();
	$min_wave=$wave_now[0];
	$max_wave=$wave_now[$nx-1];

	pgbegin(0,$dev_plot,1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.2);           # Set character height
	pgenv($min_wave,$max_wave,$y_min,$y_max,0,0);
	pgsch(0.9);           # Set character height
	pglabel("Wavelength","Flux","$title");
	pgsch(0.5);
	for ($j=0;$j<$ny;$j++) {
	    my $t=$pdl_output->slice(":,($j)");
	    pgsci(1+$j);
	    my @flux=list($t);
	    pgline($nx,\@wave_now,\@flux);
	}
	pgsci(1);
	pgclose;
	pgend;
    }
}

sub plot_test_external {
    my $pdl_wave=$_[0];
    my $pdl_flux=$_[1];
    my $pdl_model=$_[2];
    my $pdl_model_cont=$_[3];
    my ($n)=$pdl_wave->dims();
    my $i;
    open(OUT,">out_mod_res.fit_spectra");
    my $i;
    for ($i=0;$i<$n;$i++) {
	my $w=$pdl_wave->at($i);
	my $f=$pdl_flux->at($i);
	my $f_mod=$pdl_model->at($i);#*($pdl_masked->at($i));
	my $f_res=$f-$f_mod;
	my $f_cont=$pdl_model_cont->at($i);
	print OUT "$w $f $f_mod $f_red $f_cont\n";
    }
    close(OUT);
    if ($plot==1) {
	$call="plot_out_fit_mod.pl out_mod_res.fit_spectra 1/xs ";
	system($call);
    }
}


sub create_single_model() {
    my $pdl_wave=$_[0];
    my $i_now=$_[1];
    my $junk=$_[2];
    my @type=@$junk;
    my $junk=$_[3];
    my @a_c=@$junk;

    my @dims=$pdl_wave->dims();
    my $nx=$dims[0];
    my $i=0;
    my $pdl_out=zeroes($nx);
    my $j=0;

    for ($i=0;$i<$nx;$i++) {
	my $w=$pdl_wave->at($i);
	if ($type[$i_now] eq "eline") {
	    my $speed_of_light=299792.458;
	    my $factor=(1+$a_c[$i_now][3]/$speed_of_light+$a_c[$i_now][5]);
	    my $e1=1;
	    my $Y1=1;
	    if ($a_c[$i_now][2]!=0) {
		$e1=exp(-0.5*(($w-$a_c[$i_now][0]*$factor)/$a_c[$i_now][2])**2);
		$Y1=$a_c[$i_now][1]*$e1/($a_c[$i_now][2]*((2*3.1416)**0.5));
	    }
	    set($pdl_out,$i,$Y1);
	}
	if ($type[$i_now] eq "poly1d") {
	    my $ii;
	    my $Yi=0;
	    for ($ii=0;$ii<9;$ii++) {
		$Yi=$Yi+$a_c[$i_now][$ii]*($w)**($ii);
	    }
#	    print "$Yi\n";
	    set($pdl_out,$i,$Yi);
	}
    }
    return $pdl_out;
}




sub print_a_val() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_print=@$junk;
    my $junk=$_[2];
    my @e_a_print=@$junk;
    my $junk=$_[3];
    my @type=@$junk;
    my $chi_a=$_[4];
    my $ii,$j;
    print "$n_mod $chi_a\n";
    for ($ii=0;$ii<$n_mod;$ii++) {   
	print "$type[$ii] ";
	for ($j=0;$j<9;$j++) {
	    print "$a_print[$ii][$j] $e_a_print[$ii][$j] ";
	}
	print "\n";
    }
}


sub add_back_noise() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_print=@$junk;
    my $junk=$_[2];
    my @type=@$junk;
    my $chi_a=$_[3];
    my $back_noise=$_[4];
    my $ii,$j;
    for ($ii=0;$ii<$n_mod;$ii++) {
        my $e_F=$back_noise*2.354*$a_print[0][$ii][2];
        $a_print[1][$ii][1]=sqrt(($a_print[1][$ii][1])**2+$e_F**2);
    }
    return @a_print;
}




sub print_a_final() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_print=@$junk;
    my $junk=$_[2];
    my @type=@$junk;
    my $chi_a=$_[3];
    my $ii,$j;
    print "$n_mod $chi_a\n";
    for ($ii=0;$ii<$n_mod;$ii++) {   
	print "$type[$ii] ";
	for ($j=0;$j<9;$j++) {
	    print "$a_print[0][$ii][$j] $a_print[1][$ii][$j] ";
	}
	print "\n";
    }
}

sub print_a_final_file() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_print=@$junk;
    my $junk=$_[2];
    my @type=@$junk;
    my $chi_a=$_[3];
    my $outfile=$_[4];
    my $ii,$j;
    print "OUT_FILE = $outfile\n";
    open(OUT,">$outfile");
    print OUT "$n_mod $chi_a\n";
    for ($ii=0;$ii<$n_mod;$ii++) {   
	print OUT "$type[$ii] ";
	for ($j=0;$j<9;$j++) {
	    print OUT "$a_print[0][$ii][$j] $a_print[1][$ii][$j] ";
	}
	print OUT "\n";
    }
    close(OUT);
}


sub print_a_final_file_add() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_print=@$junk;
    my $junk=$_[2];
    my @type=@$junk;
    my $chi_a=$_[3];
    my $outfile=$_[4];
    my $i_val=$_[5];
    my $ii,$j;
    print "OUT_FILE = $outfile\n";
    open(OUT,">>$outfile");
#    print OUT "#ID $i_val\n";
    print OUT "$n_mod $chi_a $i_val\n";
    for ($ii=0;$ii<$n_mod;$ii++) {   
	print OUT "$type[$ii] ";
	for ($j=0;$j<9;$j++) {
	    print OUT "$a_print[0][$ii][$j] $a_print[1][$ii][$j] ";
	}
	print OUT "\n";
    }
    close(OUT);
}




#print_config_file($n_mod,$chi_sq_goal,$d_chi_sq_goal,\@a_final,\@type,\@ia,\@a0,\@a1,\@link,$out_config);

sub print_config_file() {
    my $n_mod=$_[0];
    my $chi_sq_goal=$_[1];
    my $d_chi_sq_goal=$_[2];
    my $junk=$_[3];
    my @a_print=@$junk;
    my $junk=$_[4];
    my @type=@$junk;
    my $junk=$_[5];
    my @ia=@$junk;
    my $junk=$_[6];
    my @a0=@$junk;
    my $junk=$_[7];
    my @a1=@$junk;
    my $junk=$_[8];
    my @link=@$junk;
    my $outfile=$_[9];
    my $ii,$j;
    print "OUT_CONFIG = $outfile\n";
    open(OUT,">$outfile");
    print OUT "0 $n_mod $chi_sq_goal $d_chi_sq_goal\n";
    for ($ii=0;$ii<$n_mod;$ii++) {   
	print OUT "$type[$ii]\n";
	for ($j=0;$j<9;$j++) {
	    print OUT "$a_print[0][$ii][$j] $ia[$ii][$j] $a0[$ii][$j] $a1[$ii][$j] $link[$ii][$jj]\n";
	}
    }
    close(OUT);
}


sub copy_a() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_print=@$junk;
    my @a_copy;
    my $ii,$j;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	for ($j=0;$j<9;$j++) {
	    $a_copy[$ii][$j]=$a_print[$ii][$j];
	}
    }
    return @a_copy;
}

sub copy_a_pdl() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_copy;
    my $ii,$j;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	for ($j=0;$j<9;$j++) {
	    $a_copy[$ii][$j]=$junk->at($j,$ii);
	}
    }
    return @a_copy;
}

sub random_a() {
    my $n_mod=$_[0];
    my $r_scale=$_[1];
    my $junk=$_[2];
    my @a_print=@$junk;
    my @a_copy;
    my $ii,$j;
    for ($ii=0;$ii<$n_mod;$ii++) {
	my $pdl_rnd=grandom(10);
	for ($j=0;$j<9;$j++) {	    
		my $rnd=$pdl_rnd->at($j); 
		if (($j!=3)&&($link[$ii][$j]==-1)) {		
		    $a_copy[$ii][$j]=$a_print[$ii][$j]+$r_scale*$rnd*$a_print[$ii][$j];
		} else {
		    $a_copy[$ii][$j]=$a_print[$ii][$j];
		}
	    
	}
    }
    return @a_copy;
}

#	($n_mod_fixed,\@a_fixed,\@a_type_fixed)=add_a_results($n_mod,\@a_final,\@a_type,$n_mod_fixed,\@a_fixed,\@a_type_fixed);    
sub add_a_results_elines() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_final=@$junk;
    my $junk=$_[2];
    my @a_type=@$junk;
    my $n_mod_fixed=$_[3];
    my $junk=$_[4];
    my @a_final_fixed=@$junk;
    my $junk=$_[5];
    my @a_type_fixed=@$junk;
    my $ii,$j;
    my $KK=0;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	if ($a_type[$ii] eq "eline") {
	    $a_type_fixed[$KK+$n_mod_fixed]=$a_type[$ii];	
	    for ($j=0;$j<9;$j++) {
		$a_final_fixed[$KK+$n_mod_fixed][$j]=$a_final[0][$ii][$j];
	    }
	    $KK++;
	}
    }
    $n_mod_fixed=$n_mod_fixed+$KK;
    return ($n_mod_fixed,\@a_final_fixed,\@a_type_fixed);
}

sub add_a_results() {
    my $n_mod=$_[0];
    my $junk=$_[1];
    my @a_final=@$junk;
    my $junk=$_[2];
    my @a_type=@$junk;
    my $n_mod_fixed=$_[3];
    my $junk=$_[4];
    my @a_final_fixed=@$junk;
    my $junk=$_[5];
    my @a_type_fixed=@$junk;
    my $ii,$j;
    my $KK=0;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	$a_type_fixed[$ii+$n_mod_fixed]=$a_type[$ii];	
	for ($j=0;$j<9;$j++) {
	    $a_final_fixed[$ii+$n_mod_fixed][$j]=$a_final[0][$ii][$j];
	}
    }
    $n_mod_fixed=$n_mod_fixed+$n_mod;
    return ($n_mod_fixed,\@a_final_fixed,\@a_type_fixed);
}

sub copy_a_results() {
    my $n_mod=$_[0];
    my $kk=$_[1];
    my $junk=$_[2];
    my @a_print=@$junk;
    my $junk=$_[3];
    my @a_copy=@$junk;
    my $ii,$j;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	for ($j=0;$j<9;$j++) {
	    $a_copy[$kk][$ii][$j]=$a_print[$ii][$j];
	}
    }
    return @a_copy;
}

sub mean_a_results() {
    my $n_mod=$_[0];
    my $nk=$_[1];
    my $junk=$_[2];
    my @a_print=@$junk;
    my @a_copy;
    my $ii,$j,$k;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	for ($j=0;$j<9;$j++) {
	    my @a_tmp;
	    for ($k=0;$k<$nk;$k++) { 
		$a_tmp[$k]=$a_print[$k][$ii][$j];
	    }
	    my $val=mean(@a_tmp);
	    my $e_val=sigma(@a_tmp);
	    if ($e_val==0) {
		$e_val=0.1*$val;
	    }
	    $a_copy[0][$ii][$j]=$val;
	    $a_copy[1][$ii][$j]=10*$e_val*$ia[$ii][$j];
	}
    }
    return @a_copy;
}



sub mean_a_results_last() {
    my $n_mod=$_[0];
    my $nk=$_[1];
    my $junk=$_[2];
    my @a_print=@$junk;
    my @a_copy;
    my $ii,$j,$k;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	for ($j=0;$j<9;$j++) {
	    my @a_tmp;
	    for ($k=0;$k<$nk;$k++) { 
		$a_tmp[$k]=$a_print[$k][$ii][$j];
	    }
	    my $val=mean(@a_tmp);
	    my $e_val=sigma(@a_tmp);
	    if ($e_val==0) {
		$e_val=0.1*$val;
	    }
	    $a_copy[0][$ii][$j]=$a_tmp[$nk-1];
	    $a_copy[1][$ii][$j]=$e_val*$ia[$ii][$j];
	}
    }
    return @a_copy;
}

sub mean_a_results_last_maps() {
    my $n_mod=$_[0];
    my $nk=$_[1];
    my $i_mod=$_[2];
    my $I=$_[3];
    my $J=$_[4];
    my $junk=$_[5];
    my @a_print=@$junk;
    my @a_copy;
    my $ii,$j,$k;
    for ($ii=0;$ii<$n_mod;$ii++) {   
	for ($j=0;$j<9;$j++) {
	    my @a_tmp;
	    for ($k=0;$k<$nk;$k++) { 
		$a_tmp[$k]=$a_print[$k][$ii][$j];
	    }
	    my $val=mean(@a_tmp);
	    my $e_val=sigma(@a_tmp);
	    if ($e_val==0) {
		$e_val=0.1*$val;
	    }
	    $a_copy[0][$ii][$j]=$a_tmp[$nk-1];
	    $a_copy[1][$ii][$j]=$e_val*$ia[$ii][$j];
	}
	my $II=$ii+$i_mod;
	if ($type eq "poly1d") {
	    my $sum=0;
	    my $sum_e=0;
	    my @a_wave=list($pdl_wave);
	    my $jj;
	    for ($jj=0;$jj<$#a_wave+1;$jj++) {
		for ($j=0;$j<9;$j++) {
		    $sum=$sum+$a_copy[0][$II][$j]*($a_wave[$jj])**$j;
		    $sum_e=$sum_e+$a_copy[1][$II][$j]*($a_wave[$jj])**$j;
		}
	    }
	    set($map_f,($II,$I,$J),$sum);
	    set($map_ef,($II,$I,$J),$sum_e);
	} else {	    
	    set($map_f,($II,$I,$J),$a_copy[0][$II][1]);
	    set($map_ef,($II,$I,$J),$a_copy[1][$II][1]);
	    set($map_dw,($II,$I,$J),2.354*$a_copy[0][$II][2]);
	    set($map_edw,($II,$I,$J),2.354*$a_copy[1][$II][2]);
	    set($map_v,($II,$I,$J),$a_copy[0][$II][3]);
	    set($map_ev,($II,$I,$J),$a_copy[0][$II][3]);
	}	
    }
    return @a_copy;
}

sub create_single_model_one() {
    my $pdl_wave=$_[0];
    my $i_now=$_[1];
    my $junk1=$_[2];
    my $junk=$_[3];

    my @type=@$junk1;
    my @a_c=@$junk;


    my @dims=$pdl_wave->dims();
    my $nx=$dims[0];
    my $i=0;
    my $pdl_out=zeroes($nx);
    my $j=0;

    for ($i=0;$i<$nx;$i++) {
	my $w=$pdl_wave->at($i);
	if ($type[$i_now] eq "eline") {
	    my $speed_of_light=299792.458;
	    my $factor=(1+$a_c[$i_now][3]/$speed_of_light+$a_c[$i_now][5]);
	    my $e1=1;
	    my $Y1=1;
	    if ($a_c[$i_now][2]!=0) {
		$e1=exp(-0.5*(($w-$a_c[$i_now][0]*$factor)/$a_c[$i_now][2])**2);
#		$Y1=$a_c[$i_now][1]*$e1/($a_c[$i_now][2]*((2*3.1416)**0.5));
		$Y1=1.0*$e1/($a_c[$i_now][2]*((2*3.1416)**0.5));
	    }
	    set($pdl_out,$i,$Y1);
	}
	if ($type[$i_now] eq "poly1d") {
	    my $ii;
	    my $Yi=0;
	    for ($ii=0;$ii<9;$ii++) {
		$Yi=$Yi+$a_c[$i_now][$ii]*($w)**($ii);
	    }
	    set($pdl_out,$i,$Yi);
	}
    }
    return $pdl_out;

}


sub create_single_model_poly() {
    my $pdl_wave=$_[0];
    my $ii=$_[1];
    my @dims=$pdl_wave->dims();
    my $nx=$dims[0];
    my $i=0;
    my $pdl_out=zeroes($nx);
    my $j=0;
    for ($i=0;$i<$nx;$i++) {
	my $w=$pdl_wave->at($i);
	my $Yi=($w)**($ii);
	set($pdl_out,$i,$Yi);
    }
    return $pdl_out;
}


sub smooth_ratio() {
    my $pdl_data=$_[0];
    my $pdl_model_spec_min=$_[1];
    my $sigma=$_[2];
    my $pdl_rat=$pdl_data/$pdl_model_spec_min;
    $pdl_rat->inplace->setbadtoval(1);
    if ($sigma<1) {
	$sigma=1;
    }
    my $w=int(5*2.354*$sigma);
   my @rat=list($pdl_rat);
   my @med_rat=median_filter(int(7*2.354*$sigma),\@rat);
   my $smooth_ratio=pdl(@med_rat);
    $smooth_ratio->inplace->lclip(0,1);
    $smooth_ratio->inplace->setnantobad();
    $smooth_ratio->inplace->setbadtoval(1);
    return $smooth_ratio;
}




1;




