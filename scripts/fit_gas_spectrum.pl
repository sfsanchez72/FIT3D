#!/usr/bin/perl
#
#
#




if ($#ARGV<4) {
    print "USE: fit_gas_spectrum.pl SPECTRUM.TXT LIST_OF_EMISSION_LINES min_w max_w MASK_LIST [CONFIG] [DEVICE]\n";    
    exit;
}

$spec_file=$ARGV[0];
$em_file=$ARGV[1];
$w_start=$ARGV[2];
$w_end=$ARGV[3];
$mask_list=$ARGV[4];
$no_config=0;
if ($#ARGV==5) {
    $config=$ARGV[5];
    $no_config=1;
}
if ($#ARGV==6) {
    $config=$ARGV[5];
    $device=$ARGV[6];
    $no_config=1;
}


if ($mask_list eq "none") {
    $nmask=0;
} else {
    open(FH,"<$mask_list");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$start_mask[$nmask]=$data[0];
	$end_mask[$nmask]=$data[1];
	$nmask++;
    }
    close(FH);
}

print "$nmask regions\n";


print "$spec_file $w_start $w_end\n";
$nw=0;
$y_min=10000;
$y_max=-10000;
open(FH,"<$spec_file");
open(FHOUT,">fit_spectra.input");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $id=$data[0];
    $w_test=$data[1];
    $f_test=$data[2];
    
#    print "$w_spec $w_start $w_end\n";
    if (($w_test>$w_start)&&($w_test<$w_end)) {
	$masked=0;
	for ($j=0;$j<$nmask;$j++) {
	    if (($w_test>$start_mask[$j])&&($w_test<$end_mask[$j])) {
		$masked=1;
		
	    }
	}
	if ($masked==0) {
	    print FHOUT "$id $w_test $f_test\n";
	    $w_spec[$nw]=$data[1];
	    $f_spec[$nw]=$data[2];
	    if ($y_min>$f_spec[$nw]) {
		$y_min=$f_spec[$nw];
	    }
	    if ($y_max<$f_spec[$nw]) {
	    $y_max=$f_spec[$nw];
	}
	    $nw++;
	}
    }
}
close(FHOUT);
close(FH);

if ($no_config==0) {
    open(FH,">junk.config");
    print FH "0 0 0.2 0.001\n";
    close(FH);
    #system("fit_spec_back fit_spectra.input junk.config $back_list $nback");
    if ($nf>0) {
	$call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback." ".$device;
    } else {
	$call="fit_spec_back fit_spectra.input tmp.config none 0 ".$device;
    }
    system($call);

    $n=0;
    open(FH,"<$em_file");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$w[$n]=$data[0];
	$name[$n]=$data[1];
	print "$line\n";
	$n++;
    }
    close(FH);
    print "$n lines found\n\n";
    


    print "Number of systems to include:";
    $n_s=<stdin>;
    if ($n_s!=int($n_s)) {
	print "\nError\n";
	exit;
    }
    
    if ($n_s==0) {
	print "Redshift of the system:";
	$z=<stdin>;
	chop($z);
	$vel=$z*299792.458;
	print "Systemic Velocity=$vel km/s\n";
	print "Flexibility of the velocity (km/s):";
	$flex_w=<stdin>;
	chop($flex_w);
    }    

    for ($i=0;$i<$n_s;$i++) {
	print "Redshift of the system number $i:";
	$z[$i]=<stdin>;
	chop($z[$i]);
	$vel=$z[$i]*299792.458;
	print "Systemic Velocity=$vel km/s\n";
	print "Flexibility of the velocity (km/s):";
	$flex_w=<stdin>;
	chop($flex_w);
	print "GAUSSIAN sigma of the line (Angstroms!) (GAUSS,min,max):";
	$tmp=<stdin>;
	chop($tmp);
	($dz[$i],$mdz[$i],$Mdz[$i])=split(" ",$tmp);        

    }
    
    $n_line;
    for ($i=0;$i<$n_s;$i++) {
	$nl[$i]=0;
	for ($j=0;$j<$n;$j++) {	
	    if ($j==0) {
		$fl_s[$i]=$n_line+1;
	    }
	    $w1=$w[$j]*(1+$z[$i]);
	    if (($w1>$w_start)&&($w1<$w_end)) {
		print "$name[$j] $w[$j] at $w1 found\n";
		print "Include? (y/n):";
		$include=<stdin>;
		chop($include);
		if ($include !~ "n") {
#		    $wl[$nl[$i]][$i]=$w1;
		    $wl[$nl[$i]][$i]=$w[$j];
		    print "Flux (Flux, min, Max)?:";
		    $tmp=<stdin>;
		    chop($tmp);
		    ($flux[$nl[$i]][$i],$min_flux[$nl[$i]][$i],$max_flux[$nl[$i]][$i])=split(" ",$tmp);
		    $nl[$i]++;
		    $n_line++;	
		    print "Included\n";
	    }	
	    }    
	}
	print "System $i contains $nl[$i] lines\n";
    }
    print "The total number of lines is $n_line\n";
    $n_mod=$n_line+$nback;




    $flex_dw=0.05;
    
    open(FH,">tmp.config");
    print FH "0 $n_mod 0.2 0.001\n";
    for ($i=0;$i<$n_s;$i++) {
    $j=0;
    $a=$wl[$j][$i];
    $p_a=$wl[$j][$i]-$flex_w*0.5-$dz[$i];
    $p_a2=$wl[$j][$i]+$flex_w*0.5+$dz[$i];
    $d_a=$dz[$i];
    $p_da=$mdz[$i];
    $p_da2=$Mdz[$i];
    $ff=$flux[$j][$i];
    $mff=$min_flux[$j][$i];
    $Mff=$max_flux[$j][$i];
    $vel=$z[$i]*299792.458;
    $vel1=$vel-$flex_w;
    $vel2=$vel+$flex_w;
    print FH "eline\n";
    print FH "$a\t 0\t 0\t 0\t -1\n"; 
    print FH "$ff\t 1\t $mff\t $Mff\t -1\n"; 
    print FH "$d_a\t 1\t $p_da\t $p_da2\t -1\n"; 
    print FH "$vel\t 1\t $vel1\t $vel2\t -1\n"; 
    for ($k=0;$k<5;$k++) {
	print FH "0\t 0\t 0\t 0\t -1\n";
    }
    for ($j=1;$j<$nl[$i];$j++) {
	$a=$wl[$j][$i];
	$p_a=$wl[$j][$i]-$wl[0][$i];
	$p_a2=0.0;
	$d_a=$dz[$i];
	$p_da=$mdz[$i];
	$p_da2=$Mdz[$i];
	$ff=$flux[$j][$i];
	$mff=$min_flux[$j][$i];
	$Mff=$max_flux[$j][$i];
	print FH "eline\n";
	print FH "$a\t 0\t 0\t 0\t -1\n"; 
	print FH "$ff\t 1\t $mff\t $Mff\t -1\n"; 
#	print FH "$d_a\t 1\t $p_da\t $p_da2\t 1\n"; 
	print FH "$d_a\t 1\t 0\t 0\t $fl_s[$i]\n"; 
	print FH "$vel\t 1\t 0\t 0\t $fl_s[$i]\n"; 
	for ($k=0;$k<5;$k++) {
	    print FH "0\t 0\t 0\t 0\t -1\n";
	}
    }


}

    for ($j=0;$j<$nback;$j++) {
	print "Level of the Background $j/$nback (level, min,max)?";
	$tmp=<stdin>;
	chop($tmp);
	if ($tmp!="") {
	    ($level,$lmin,$lmax)=split(" ",$tmp);
	} else {
	    print "$level,$lmin,$lmax\n";
	}
	print "GAUSSIAN sigma of the line (Angstroms!) (GAUSS,min,max):";
	$tmp=<stdin>;
	chop($tmp);
	if ($tmp!="") {
	    ($dz,$mdz,$Mdz)=split(" ",$tmp);        
	} else {
	    print "$dz,$mdz,$Mdz\n";
	}
	print FH "back\n";
	print FH "$j\t 0\t 0\t 0\t -1\n"; 
	print FH "$level\t 1\t $lmin\t $lmax\t -1\n";
	if (($n_s==0)&&($j==0)) {
	    $vel1=$vel-$flex_w;
	    $vel2=$vel+$flex_w;
	    print FH "$vel\t 1\t $vel1\t $vel2\t -1\n"; 
	    print FH "$dz\t 1\t $mdz\t $Mdz\t -1\n";
	} else {
	    print FH "$vel\t 1\t 0\t 0\t 1\n"; 
	    print FH "$dz\t 1\t 0\t 0\t 1\n";

	}
	for ($k=0;$k<5;$k++) {
	    print FH "0\t 0\t 0\t 0\t -1\n";
	}
    }
close(FH);

} else {
    system("cp $config tmp.config");
}

#system("rm out.model");
#system("rm model.fits");
#system("rm res.fits");
system("rm out_config.fit_spectra");
system("rm out.fit_spectra");
if ($#ARGV==3) {
 system("emacs tmp.config &");
 print "Make the changes by hand in the Config file\n";
 print "Save it and press <Enter>\n";
 <stdin>;
}
    if ($nf>0) {
	$call="fit_spec_back_fits fit_spectra.input tmp.config BACK_JUNK.fits ".$nback." ".$device;
    } else {
	$call="fit_spec_back fit_spectra.input tmp.config none 0 ".$device;
    }
    system($call);

#$call="fit_spec_back fit_spectra.input tmp.config ".$back_list." ".$nback;
#system($call);
print "$call\n";
print "Results saved in 'out.fit_spectra' and 'out_config.fit_spectra'\n";
print "For re-running the program, use:\n fit_spec_back fit_spectra.input out_config.fit_spectra $back_list $nback\n"; 
#print "For re-running the program, use:\n fit_spectra fit_spectra.input tmp.config\n"; 

$n=0;
open(FH,"<out_mod_res.fit_spectra");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $cont[$n]=$data[3];
    $n++;
}
close(FH);
$cont_level=mean(@cont);
$s_cont_level=sigma(@cont);
print "CONTINUUM=$cont_level+-$s_cont_level\n";
#
# Once done the analysis, we need to create
# the position Intensity, Velocity and Dispersion Maps.
#




exit;

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

sub sigma {
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

  return $stddev;
}

