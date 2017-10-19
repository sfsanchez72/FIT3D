#!/usr/bin/perl
#
#
#
#use PGPLOT;


# This program $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #s the and output file
# like the one created by "fit_3D_fg.pl"
# and fits over the datacube, using the
# template file (CONFIG_FILE), in the 
#
#
#


$tk_e3d="/home/sanchez/sda2/code/lyon/v3d/user/bin/tk_e3d";
if ($#ARGV<3) {
    print "USE: fit_3D_force.pl E3D_FITSFILE CONFIG_FILE INPUT_FILE OUT_FILE [FRE_FRACT] [BIAS]\n";    
    exit;
}

 $e3d_file=$ARGV[0];
 $config_file=$ARGV[1];
 $config_tmp="config_tmp.fit_3D_force";
 $input_file=$ARGV[2];
 $out_file=$ARGV[3];
    $free_frac_def=0;
$bias=0;
if ($#ARGV>3) {
    $free_frac=$ARGV[4];
    $free_frac_def=1;
    $bias=$ARGV[5];
}
$auto="n";
$slice_file="slice_fit_3D_fg.tmp";


#
# We read the config file
#
#

print "Reading the initical config file: $input_file\n";
open(FH,"<$config_file");
$header=<FH>;
chop($header);
($junk,$n_mod,$chi_sq,$dchi_sq)=split(" ",$header);
srand;
for ($i=0;$i<$n_mod;$i++) {
    $line=<FH>; chop($line);
    $model[$i]=$line;
    for ($j=0;$j<9;$j++) {
	$line=<FH>; chop($line);
	@data=split(" ",$line);
	$a[$i][$j]=$data[0];
	$ia[$i][$j]=$data[1];
	$a_m[$i][$j]=$data[2];
	$a_M[$i][$j]=$data[3];
	$flag[$i][$j]=$data[4];
    }
}
close(FH);
print "DONE\n";




#system("emacs $config_file &");	

#$vis=$ARGV[4];

#system("rm $out_file");
if (($vis!=0)&&($vis!=1)) {
    $vis=1;
}


print "Reading the data from the image\n";
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
print "DONE\n";

#print "$npix, $nb_spec, $start_w\n";
#exit;

system("rm org.fit_3D_fg");
system("rm mod.fit_3D_fg");
system("rm res.fit_3D_fg");

#
# Loop over all the spectra!!!
#

open(INPUT,"<$input_file");
for ($js=1;$js<=$npix;$js++) {
#
# We extract the slice number $js
#
    $j_end=$js+1;

#
# Reads the data
#
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
# We create the proper $config_tmp, using the entries on the
# INPUT file and the config file
#
    $input=<INPUT>; chop($input);
    if ($input!=$n_mod) {
	print "******************** ERROR *******************************\n";
	print "The number of models in the file '$input_file' is different\n";
	print "to the number of models in the file '$config_file'\n";	
	exit;
    }
    for ($i=0;$i<$n_mod;$i++) {
	$input=<INPUT>; chop($input);
	@data=split(" ",$input);
	for ($j=0;$j<9;$j++) {
	    $a[$i][$j]=$data[1+2*$j];
	    if (($ia[$i][$j]==1)&&($flag[$i][$j]==-1)&&($free_frac_def==1)) {
		$a_m[$i][$j]=$a[$i][$j]-$free_frac*abs($a[$i][$j]);
		$a_M[$i][$j]=$a[$i][$j]+$free_frac*abs($a[$i][$j]);
	    }
#    print "**** $a[$i][$j] ";
	    if ($j==2) {
#		$a[$i][$j]=$data[1+2*$j]+((-1)**(int(0.5+rand(1))))*rand(0.5)*$data[1+2*$j];
	    }
#    print "$a[$i][$j] $i $j\n";
#	    print "$a[$i][$j] ";
	}
#	print "\n";
#	print "$input\n";
    }

    open(TMP,">$config_tmp");
    print TMP "$junk $n_mod $chi_sq $dchi_sq\n";
    for ($i=0;$i<$n_mod;$i++) {
	print TMP "$model[$i]\n";
#	print "IN =$model[$i]\n";
	for ($j=0;$j<9;$j++) {
	    print TMP "$a[$i][$j] $ia[$i][$j] $a_m[$i][$j] $a_M[$i][$j] $flag[$i][$j]\n";
#	    print "IN =$a[$i][$j] $ia[$i][$j] $a_m[$i][$j] $a_M[$i][$j] $flag[$i][$j]\n";
	}
    }
    close(TMP);
    
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
#	print "OUT=$line";
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


#    if ($wait==1) {
#	print "PRESS a KEY\n";
#	<stdin>;
#    }
}
close(INPUT);


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

