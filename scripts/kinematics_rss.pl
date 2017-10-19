#!/usr/bin/perl
#
#
#
#use PGPLOT;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");
$tk_e3d="/home/sanchez/sda2/code/lyon/v3d/user/bin/tk_e3d";
if ($#ARGV<5) {
    print "USE: kinematics.pl RSS_FITSFILE LINE_FILE START_WAVELENGTH END_WAVELENGTH OUT_FILE DEVICE\n";    
    exit;
}

$rss_file=$ARGV[0];
$em_file=$ARGV[1];
$w_start=$ARGV[2];
$w_end=$ARGV[3];
$spec_file="spectrum.kin";
$slice_file="slice.kin";
$out_file=$ARGV[4];
$device=$ARGV[5];
$auto="n";
#$vis=$ARGV[4];

#system("rm $out_file");
if (($vis!=0)&&($vis!=1)) {
    $vis=1;
}


#open(TCL,">kin.tcl");
#print TCL "create_env euro3d /null\n";
#print TCL "euro3d load_file $e3d_file\n";
#print TCL "set data [euro3d ask_e3d_info]\n";
#print TCL "set Fid [open junk.txt {RDWR CREAT} 0777]\n";
#print TCL "puts \$Fid \"\$data\"\n";
#print TCL "euro3d create_slice 0 1 1 0\n";
#print TCL "euro3d save_slice 1 $slice_file Flux 1\n";
#print TCL "close \$Fid\n";
#print TCL "exit\n";
#close(TCL);
#system("$tk_e3d kin.tcl");




#open(OUT,"<junk.txt");
#$line=<OUT>;
#chop($line);
#@data=split(" ",$line);

#
# We read the RSS image
#
my $nt=read_naxes($rss_file);
($npix,$nb_spec) = @$nt;
my @rss=read_img($rss_file);
($start_w,$delta_w)=read_img_headers($rss_file,["CRVAL1","CDELT1"]);

print "$nb_spec spectra at $rss_file\n";
#print "$npix, $nb_spec, $start_w\n";
#exit;


#
# Loop over all the spectra!!!
#
for ($js=0;$js<$nb_spec;$js++) {
    $nw=0;
    $y_min=10000;
    $y_max=-10000;
#    open(FH,"<$spec_file");
    open(FHOUT,">fit_spectra.input");
    for ($i=0;$i<$npix;$i++) {	
	$id=$k+1;
	$w_test=$start_w+$delta_w*$i;
	$f_test=$rss[$js][$i];
#	print "$w_test $f_test\n";
	if (($w_test>$w_start)&&($w_test<$w_end)) {
	    print FHOUT "$id $w_test $f_test\n";
	    $w_spec[$nw]=$w_test;
	    $f_spec[$nw]=$f_test;
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
#    close(FH);
    $w_delta_spec=$w_spec[1]-$w_spec[0];

#
# Ask only for the first spectrum
#
    if ($js==0) {
	
	open(FH,">junk.config");
	print FH "0 0 0.2 0.001\n";
	close(FH);
#	system("fit_spectra fit_spectra.input junk.config");
	system("fit_spectra fit_spectra.input junk.config $device");
	
	
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
	$back_pa=0.1*$back;
	$back_pa2=10**$back;
	
	print FH "poly1d\n";
	print FH "$back\t 1\t $back_pa\t $back_pa2\t -1\n"; 
	for ($k=0;$k<$ord_back;$k++) {
	    printf FH "0.0001\t 1\t -5.0\t 5.0\t -1\n";
	}
	for ($k=$ord_back;$k<(8-$ord_back-1);$k++) {
	    printf FH "0\t 0\t 0\t 0\t -1\n";
	}
	
	close(FH);

	system("rm out_config.fit_spectra");
	system("rm out.fit_spectra");
	system("emacs tmp.config &");	
	print "Make the changes by hand in the Config file\n";
	print "Save it and press <Enter>\n";
	<stdin>;
    }

    if ($auto !~ "y") {
	$command="r";
	$n_fit=1;
	while ($command !~ "s") {
	    if ($command =~ "r") {
		print "Fitting Num. $n_fit of the spectrum $js (of $nb_spec)\n";
		$call="fit_spectra fit_spectra.input tmp.config ".$device;
		system($call);
		$n_fit++;
	    }
	    
	    print "Options:\n";
	    print "[s] save the results\n";
	    print "[r] repeat the fitting\n";
	    print "[c] estimate continuum\n";
	    print "[m] copy the output of the last fit to the new config\n";
	    print "[a] Set the automatic fitting\n\n";
	    $command=<stdin>;
	    chop($command);
	    
	    if ($command =~ "a") {
		$command="s";
		$auto="y";
	    }

	    if ($command =~ "m") {
		system("cp out_config.fit_spectra tmp.config");
		system("emacs tmp.config &");
	    }
	    
	    if ($command =~ "c") {
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
		
	    }
	    
	}
    } else {
	print "Fitting spectrum $js (of $nb_spec)\n";
	$call="fit_spectra fit_spectra.input tmp.config ".$device;
	system($call);
    }
#
# We copy the values derived from the fitting
#
    open(FH,"<out.fit_spectra");
    if ($js==0) {
	open(FHOUT,">$out_file");
    } else {
	open(FHOUT,">>$out_file");
    }
    while($line=<FH>) {
	print FHOUT "$line";
    }
    close(FHOUT);
    close(FH);
#
# We copy the model and the residual to an ascci file
#
    open(FH,"out_mod_res.fit_spectra");
    if ($js==0) {
	open(ORG,">org.kinematics");
	open(MOD,">mod.kinematics");
	open(RES,">res.kinematics");
	print ORG "# $nw $nb_spec $w_spec[0] $w_delta_spec\n";
	print MOD "# $nw $nb_spec $w_spec[0] $w_delta_spec\n";
	print RES "# $nw $nb_spec $w_spec[0] $w_delta_spec\n";
    } else {
	open(ORG,">>org.kinematics");
	open(MOD,">>mod.kinematics");
	open(RES,">>res.kinematics");
    }
    $n_line=0;
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	if ($n_line==0) {
	    $start_out_w=$data[0];
	}
	if ($n_line==1) {
	    $delta_out_w=$data[0];
	}
	print ORG "$data[1]\n";
	print MOD "$data[2]\n";
	print RES "$data[3]\n";
	$n_line++;
    }
    close(ORG);
    close(MOD);
    close(RES);
    close(FH);
}

print "$result saved at $out_file\n";

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

