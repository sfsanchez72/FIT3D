#!/usr/bin/perl
#
#
#
#use PGPLOT;



#$tk_e3d="/home/sanchez/sda2/code/lyon/v3d/user/bin/tk_e3d";
$tk_e3d="tk_e3d";
if ($#ARGV<3) {
    print "USE: fit_3D_fg.pl E3D_FITSFILE CONFIG_FILE OUT_FILE [BIAS]\n";    
    exit;
}

$e3d_file=$ARGV[0];
$config_file=$ARGV[1];
$config_tmp=$config_file;
$out_file=$ARGV[2];
$bias=$ARGV[3];
$auto="n";
$slice_file="slice_fit_3D_fg.tmp";

#system("emacs $config_file &");	

#$vis=$ARGV[4];

#system("rm $out_file");
if (($vis!=0)&&($vis!=1)) {
    $vis=1;
}


open(TCL,">fit.tcl");
print TCL "create_env euro3d /null\n";
print TCL "euro3d load_file $e3d_file\n";
print TCL "set data [euro3d ask_e3d_info]\n";
print TCL "set Fid [open junk.txt {RDWR CREAT} 0777]\n";
print TCL "puts \$Fid \"\$data\"\n";
print TCL "euro3d create_slice 0 1 1 0\n";
print TCL "euro3d save_slice 1 $slice_file Flux 1\n";
print TCL "close \$Fid\n";
print TCL "exit\n";
close(TCL);

system("$tk_e3d fit.tcl");


open(OUT,"<junk.txt");
$line=<OUT>;
chop($line);
@data=split(" ",$line);
$npix=$data[0];
$nb_spec=$data[1];
$start_w=$data[2];
$end_w=$data[3];
$delta_w=$data[4];
close(OUT);

print "$npix channels ($start_w-$end_w)  at $e3d_file\n";
#print "$npix, $nb_spec, $start_w\n";
#exit;

system("rm org.fit_3D_fg");
system("rm mod.fit_3D_fg");
system("rm res.fit_3D_fg");

#
# Loop over all the spectra!!!
#
for ($js=1;$js<=$npix;$js++) {
#
# We extract the slice number $js
#
    $j_end=$js+1;

    system("rm $slice_file\n");
    open(TCL,">fit.tcl");
    print TCL "create_env euro3d /null\n";
    print TCL "euro3d load_file $e3d_file\n"; 
    print TCL "euro3d create_slice $js $j_end 1 0\n";
    print TCL "euro3d save_slice 1 $slice_file Flux 1\n";
    print TCL "exit\n";
    close(TCL);
    system("$tk_e3d fit.tcl");
    


    system("rm out_config.fit_2D_map");
    system("rm out.fit_2D_map");

#
# The first guess 
# we use an unique configurantion file
#
    $call="fit_2D_map ".$slice_file." ".$config_tmp." ".$bias." |";
    open(CALL,$call);
    $fit=0;
    while($line_out=<CALL>) {	
	if ($line_out =~ "Singular") {
	    $fit=1;
	}
	print "$line_out";
    }
    close(CALL);
#    system($call);

#
# Force
#
#    $fit=0;
    if ($fit==0) {
	$call="cp out_config.fit_2D_map config.fit_2D.tmp";
	system($call);
	$config_tmp="config.fit_2D.tmp";
	$call="cp out_config.fit_2D_map last_cf_fit_2D.tmp";
	system($call);
    } else {
#$config_file=$ARGV[1];
	$call="cp ".$config_file." config.fit_2D.tmp";
	system($call);
	$config_tmp="config.fit_2D.tmp";
    }

#
# We copy the values derived from the fitting
#
    open(FH,"<out.fit_2D_map");
    if ($js==1) {
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
    if ($js==1) {
	open(ORG,">org.fit_3D_fg");
	open(MOD,">mod.fit_3D_fg");
	open(RES,">res.fit_3D_fg");
	print ORG "# $npix $nb_spec $start_w $end_w $delta_w\n";
	print MOD "# $npix $nb_spec $start_w $end_w $delta_w\n";
	print RES "# $npix $nb_spec $start_w $end_w $delta_w\n";
    } else {
	open(ORG,">>org.fit_3D_fg");
	open(MOD,">>mod.fit_3D_fg");
	open(RES,">>res.fit_3D_fg");
    }
    $n_line=0;


    open(ORG2D,"<$slice_file");
    while($line=<ORG2D>) {
	chop($line);
	@data=split(" ",$line);
	print ORG "$data[3]\n";
	$n_line++;
    }
    close(MOD2D);

    open(MOD2D,"<fit_2D_map.mod");
    while($line=<MOD2D>) {
	chop($line);
	@data=split(" ",$line);
	print MOD "$data[3]\n";
	$n_line++;
    }
    close(MOD2D);

    open(RES2D,"<fit_2D_map.res");
    while($line=<RES2D>) {
	chop($line);
	@data=split(" ",$line);
	print RES "$data[3]\n";
	$n_line++;
    }
    close(RES2D);

    close(MOD);
    close(RES);


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

