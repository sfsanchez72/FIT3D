#!/usr/bin/perl
#
#
#

#use PGPLOT;

if ($#ARGV<3) {
    print "USE:fit_spectra.pl SPECTRUM.TXT LINE_FILE START_WAVELENGTH END_WAVELENGTH [CONFIG]\n";    
    exit;
}

$spec_file=$ARGV[0];
$em_file=$ARGV[1];
$w_start=$ARGV[2];
$w_end=$ARGV[3];
if ($#ARGV==4) {
    $config=$ARGV[4];
}


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
close(FHOUT);
close(FH);

if ($#ARGV==3) {

    open(FH,">junk.config");
    print FH "0 0 0.2 0.001\n";
    close(FH);
    system("fit_spectra fit_spectra.input junk.config");


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
    
    for ($i=0;$i<$n_s;$i++) {
	print "Redshift of the first system:";
	$z[$i]=<stdin>;
	chop($z[$i]);
	print "Flexibility on the Redshift (Angstroms):";
	$flex_w=<stdin>;
	chop($flex_w);
	print "Velocity dispersion (Angstroms!) (vel,min,max):";
	$tmp=<stdin>;
	chop($tmp);
	($dz[$i],$mdz[$i],$Mdz[$i])=split(" ",$tmp);        
	print "Flux level of the first system:";
	$f[$i]=<stdin>;
	chop($f[$i]);
	print "MIN Flux level of the first system:";
	$f_m[$i]=<stdin>;
	chop($f_m[$i]);
	print "MAX Flux level of the first system:";
	$f_M[$i]=<stdin>;
	chop($f_M[$i]);

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
		    $wl[$nl[$i]][$i]=$w1;
		    print "Flux (Flux, min, Max)? (y/n):";
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
    $n_mod=$n_line+1;




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
    print FH "gauss1d\n";
    print FH "$a\t 1\t $p_a\t $p_a2\t -1\n"; 
    print FH "$ff\t 1\t $mff\t $Mff\t -1\n"; 
    print FH "$d_a\t 1\t $p_da\t $p_da2\t -1\n"; 
    for ($k=0;$k<6;$k++) {
	printf FH "0\t 0\t 0\t 0\t -1\n";
    }
    for ($j=1;$j<$nl[$i];$j++) {
	$a=$wl[0][$i];
	$p_a=$wl[$j][$i]-$wl[0][$i];
	$p_a2=0.0;
	$d_a=$dz[$i];
	$p_da=$mdz[$i];
	$p_da2=$Mdz[$i];
	$ff=$flux[$j][$i];
	$mff=$min_flux[$j][$i];
	$Mff=$max_flux[$j][$i];
	print FH "gauss1d\n";
	print FH "$a\t 1\t $p_a\t $p_a2\t $fl_s[$i]\n"; 
	print FH "$ff\t 1\t $mff\t $Mff\t -1\n"; 
#	print FH "$d_a\t 1\t $p_da\t $p_da2\t 1\n"; 
	print FH "$d_a\t 1\t 0\t 0\t $fl_s[$i]\n"; 
	for ($k=0;$k<6;$k++) {
	    printf FH "0\t 0\t 0\t 0\t -1\n";
	}
    }
}
    
    print "BACKGROUND order:";
    $ord_back=<stdin>;
    chop($ord_back);
    if ($ord_back>8) { 
	$ord_back=8;
    }
    if ($ord_back<0) {
	$ord_back=0;
    }
    
print "BACKGROUND LEVEL:";
$back=<stdin>;
chop($back);
$back_pa=-100*$back;
$back_pa2=100**$back;

print FH "poly1d\n";
print FH "$back\t 1\t $back_pa\t $back_pa2\t -1\n"; 
for ($k=0;$k<$ord_back;$k++) {
    printf FH "0.0001\t 1\t -5.0\t 5.0\t -1\n";
}
for ($k=$ord_back;$k<(8-$ord_back-1);$k++) {
    printf FH "0\t 0\t 0\t 0\t -1\n";
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
$call="fit_spectra fit_spectra.input tmp.config";
system($call);
print "Results saved in 'out.fit_spectra' and 'out_config.fit_spectra'\n";
print "For re-running the program, use:\n fit_spectra fit_spectra.input out_config.fit_spectra\n"; 
print "For re-running the program, use:\n fit_spectra fit_spectra.input tmp.config\n"; 

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

