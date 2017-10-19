#!/usr/bin/perl
#
#
#
#use PGPLOT;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda1/perl/MY/my.pl");
#$tk_e3d="/home/sanchez/sda2/code/lyon/v3d/user/bin/tk_e3d";
$tk_e3d="tk_e3d";
if ($#ARGV<4) {
    print "USE: fit_3D_fg.pl RSS_FILE PT_FILE CONFIG_FILE OUT_FILE BIAS [Threshold]\n";    
    exit;
}

$rss_file=$ARGV[0];
$pt_file=$ARGV[1];
$config_file=$ARGV[2];
$config_tmp=$config_file;
$out_file=$ARGV[3];
$bias=$ARGV[4];
$th=-1e100;
if ($#ARGV==5) {
    $th=$ARGV[5];
}
$auto="n";
$slice_file="slice_fit_3D_fg.tmp";

#
# We read the PT
#
$ns=0;
open(FH,"<$pt_file");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    if (($line!="R")&&($line!="C")&&($line!="H")&&($line!="S")) {
	$id[$ns]=$data[0];
	$x[$ns]=$data[1];
	$y[$ns]=$data[2];
	$ns++;
    }
}
close(FH);


#system("emacs $config_file &");	

#$vis=$ARGV[4];

#system("rm $out_file");
if (($vis!=0)&&($vis!=1)) {
    $vis=1;
}

my $nt=read_naxes($rss_file);
($npix,$nb_spec) = @$nt;
my @rss=read_img($rss_file);
($start_w,$delta_w)=read_img_headers($rss_file,["CRVAL1","CDELT1"]);

print "$npix channels ($start_w-$end_w)  at $rss_file\n";
if ($ns!=$nb_spec) {
    print "$ns!=$nb_spec (Error)\n";
    exit;
}

#print "Press enter to start to fit <ENTER>\n";

#<stdin>;

system("rm org.fit_3D_fg");
system("rm mod.fit_3D_fg");
system("rm res.fit_3D_fg");

#
# Loop over all the spectral pixels!!!
#
for ($js=0;$js<$npix;$js++) {
#
# We extract the slice number $js
#
    $j_end=$js+1;

#    system("rm $slice_file\n");
#    open(TCL,">fit.tcl");
#    print TCL "create_env euro3d /null\n";
#    print TCL "euro3d load_file $rss_file\n"; 
#    print TCL "euro3d create_slice $js $j_end 1 0\n";
#    print TCL "euro3d save_slice 1 $slice_file Flux 1\n";
#    print TCL "exit\n";
#    close(TCL);
#    system("$tk_e3d fit.tcl");
 
    open (FH,">$slice_file");
    for ($j=0;$j<$nb_spec;$j++) {
	$k=$j+1;
	if ($rss[$j][$js]>$th) {
	    print FH "$k $x[$j] $y[$j] $rss[$j][$js]\n";
	}
    }
    close(FH);


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
    $k=0;
    while($line=<ORG2D>) {
	chop($line);
	@data=split(" ",$line);
	$a_org[$k][$js]=$data[3];
	print ORG "$data[3]\n";
	$k++;
	$n_line++;
    }
    close(ORG2D);

    $k=0;
    open(MOD2D,"<fit_2D_map.mod");
    while($line=<MOD2D>) {
	chop($line);
	@data=split(" ",$line);
	print MOD "$data[3]\n";
	$a_mod[$k][$js]=$data[3];
	$k++;
	$n_line++;
    }
    close(MOD2D);

    $k=0;
    open(RES2D,"<fit_2D_map.res");
    while($line=<RES2D>) {
	chop($line);
	@data=split(" ",$line);
	$a_res[$k][$js]=$data[3];
	print RES "$data[3]\n";
	$k++;
	$n_line++;
    }
    close(RES2D);

    close(MOD);
    close(RES);


}

print "$result saved at $out_file\n";

#
# Saving the fits files...
#
system("rm mod_fit_3D_fg.fits");
write_rss("mod_fit_3D_fg.fits",$npix,$nb_spec,$start_w,$delta_w,\@a_mod);
system("rm res_fit_3D_fg.fits");
write_rss("res_fit_3D_fg.fits",$npix,$nb_spec,$start_w,$delta_w,\@a_res);



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

