#!/usr/bin/perl
#
use Statistics::OLS;
use Math::FFT;
use Math::Stat;
use Math::Spline qw(spline linsearch binsearch);
use Math::Derivative qw(Derivative2);
use Math::Approx;
use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;
use PDL::Slatec;
use PDL::Image2D;
#use PDL::Matrix;


$F3D=$ENV{FIT3D_PATH}; $F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("$F3D/scripts/my.pl");


if ($#ARGV<0) {
    print "USE: auto_redu.pl DATE(YYYYMMDD) GLUE[0=NO/1=YES] [DX DY spots] [ARCHIVE CALIB]\n";
    exit;
}
$TODAY=$ARGV[0];
$GLUE=$ARGV[1];
$DX=$ARGV[2];
$DY=$ARGV[3];
$CALIB=0;
if ($#ARGV==4) {
    $CALIB=$ARGV[4];
}




$TOMORROW=$TODAY+1;

$DIR_CONF="/disk-b/sanchez/ppak/legacy/";

open(LOG,">auto_redu_arc_std.log");
open(PROD,">auto_redu_arc_std.out");




#
# We glue all the CCDs
#
if ($GLUE==1) {
    print "Glueing the CCD data...\n";
    open(DIR,"ls *a.fits |");
    while($file=<DIR>) {
	chop($file);
	$cut="a.fits";
	$file =~ s/$cut//;
	$call="glue_new_pmas.pl ".$file;
	system($call);
    }
    close(DIR);
    print "DONE\n";
}

$call="cp ".$DIR_CONF."/CCDFlat.fits .";
mycall($call);
$call="cp ".$DIR_CONF."/ARC_V1200.dist.id ARC.dist.id";
mycall($call);
$call="cp ".$DIR_CONF."/ARC_V1200.disp.id ARC.disp.id";
mycall($call);
$call="cat /meteo/CAVEX/cavex.dat | grep ".$TODAY." > cavex.dat";
mycall($call);
$call="cat /meteo/CAVEX/cavex.dat | grep ".$TOMORROW." >> cavex.dat";
mycall($call);

#exit;

$nc=0;
open(FH,"<cavex.dat");
while($line=<FH>) {
    chop($line);
    @data=split(" ",$line);
    $date[$nc]=$data[1].$data[2];
    $av[$nc]=$data[5];
    $nc++;
}
close(FH);

#exit;



if ($GLUE==1) {

$call="get_object_V600.pl > object.txt";
mycall($call);



#
# We measure the spots location
#
$call="spots_list.pl object.txt ".$DIR_CONF."/spots_V1200.short ".$DX." ".$DY;
mycall($call);
$call="plot_spots_all.pl spots_".$TODAY.".ps/CPS 200";
mycall($call);



$call="get_object_ext.pl > object_ext.txt";
mycall($call);

}




#
# How many targets are?
#
$i=1;
$n_unc=0;
$n_com=0;
do {
    $OBJ_I="obj_".$i;
    open(FH,"cat object.txt | grep ' $OBJ_I' | wc -l |");
    $nfound=<FH>;
    chop($nfound);
    close(FH);
#    if (($nfound<9)&&($nfound!=0)) {
    if ($nfound>0) {
	if ($nfound<9) {
	    $Jobj[$n_unc]=$i;
	    $n_unc++;
	} else {
	    $n_com++;
	}
    }
#    print "$FFF\n";
#    print "$i $nfound $n_com $n_unc\n";
    $i++;
} while ($nfound!=0);






#print "Nunc=$n_unc Ncom=$n_com\n";

#exit;



$call="cat object.txt | grep BIAS > bias.txt";
mycall($call);

open(BIAS,"<bias.txt");
open(OUT,">bias.out");
while($line=<BIAS>) {
    chop($line);
    @data=split(" ",$line);
    $file=$data[0];
    open(IMSTAT,"imstats.pl $file 520 1080 50 50 | ");
    while($imstat=<IMSTAT>) {
	chop($imstat);
	@dat_imstat=split(" ",$imstat);
	$mean=$dat_imstat[0];
	$rms=$dat_imstat[1];
	if ((abs($mean-653)<4)&&($rms<3.7)) {
	    print OUT "$file\n";
	}
    }
    close(IMSTAT);
}
close(OUT);
close(BIAS);

$call="wc -l bias.out > nbias.txt";
mycall($call);

open(BIAS,"<nbias.txt");
$line=<BIAS>;
chop($line);
@data=split(" ",$line);
$nbias=$data[0];
close(BIAS);



#exit;


#
# Now we reduce the calibration stars
#
$call="cat object.txt | grep ' std_1' > std_1.txt";
mycall($call);
$nstd_1=0;
open(STD_1,"<std_1.txt");
while($line=<STD_1>) {
    chop($line);
    @data=split(" ",$line);
    $std_1=$data[0];
    $name_std_1=$data[1];
    $cut="std_1_";
    $name_std_1 =~ s/$cut//;
    $nstd_1++;
}
close(STD_1);


if ($CALIB==0) {

if ($nstd_1>0) {

    $call="cat object.txt | grep ' cont_std_1' > cont_std_1.txt";
    mycall($call);
    open(CONT_STD_1,"<cont_std_1.txt");
    while($line=<CONT_STD_1>) {
	chop($line);
	@data=split(" ",$line);
	$cont_std_1=$data[0];
    }
    close(CONT_STD_1);
    
    $call="cat object.txt | grep ' arc_std_1' > arc_std_1.txt";
    mycall($call);
    open(ARC_STD_1,"<arc_std_1.txt");
    while($line=<ARC_STD_1>) {
	chop($line);
	@data=split(" ",$line);
	$arc_std_1=$data[0];
    }
    close(ARC_STD_1);
    
    $call="reduce_all.pl ".$std_1." ".$cont_std_1." ".$arc_std_1."  1 ".$DIR_CONF."/reduce.flat.V1200.config > reduce_junk.log";
    mycall($call);

    $pre_std_1=$std_1;
    $pre_std_1 =~ s/.fits//;
#($pre_std_1,$post)=split(/\./,$std_1);
    
    $call="radial_sum_rss.pl ".$pre_std_1.".sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 std_1.rss.fits";
    mycall($call);
    $call="img2spec.pl std_1.rss.fits 1 spec_std_1.txt";
    mycall($call);
    
    open(EXPTIME,"read_img_header.pl $std_1 EXPTIME,AIRMASS,DATE-OBS |");
    $line=<EXPTIME>;
    ($file,$exptime,$airmass,$date_obs)=split(" ",$line);
    close(EXPTIME);
#print "$name_std_1 $exptime $airmass\n";
    $factor=1/($exptime); # Filling-factor!
    $date_obs =~ s/T//;
    $date_obs =~ s/\-//g;
    $date_obs =~ s/\://g;
    $extinction=0.15;
    $NC=0;
	while (($date[$NC]<$date_obs)&&($NC<$nc)) {
#    while ($date[$NC]<$date_obs) {
	if ($av[$NC]>$extinction) {
	    $extinction=$av[$NC];
	}
    $NC++;
    }


    $call="calib_det.pl spec_std_1.txt ".$DIR_CONF."/f".$name_std_1.".dat ratio_std1.txt ".$factor." ".$extinction." ".$airmass." 31 15 0 3800 5000";
    mycall($call);


#
# We check if there is a second calibration start, to test the accuracy of the calibration
#

    $n_std2=0;
    $call="cat object.txt | grep ' std_2' > std_2.txt";
    mycall($call);
    open(STD_2,"<std_2.txt");
    while($line=<STD_2>) {
	chop($line);
	@data=split(" ",$line);
	$std_2=$data[0];
	$name_std_2=$data[1];
	$cut="std_2_";
	$name_std_2 =~ s/$cut//;
	$n_std2++;
    }
    close(STD_2);
    
    if ($n_std2>0) {
	$call="cat object.txt | grep ' cont_std_2' > cont_std_2.txt";
	mycall($call);
	open(CONT_STD_2,"<cont_std_2.txt");
	while($line=<CONT_STD_2>) {
	    chop($line);
	    @data=split(" ",$line);
	    $cont_std_2=$data[0];
	}
	close(CONT_STD_2);
	
	$call="cat object.txt | grep ' arc_std_2' > arc_std_2.txt";
	mycall($call);
	open(ARC_STD_2,"<arc_std_2.txt");
	while($line=<ARC_STD_2>) {
	    chop($line);
	    @data=split(" ",$line);
	    $arc_std_2=$data[0];
	}
	close(ARC_STD_2);

    $call="reduce_all.pl ".$std_2." ".$cont_std_2." ".$arc_std_2."  1 ".$DIR_CONF."/reduce.flat.V1200.config > reduce_junk.log";
    mycall($call);

    $pre_std_2=$std_2;
    $pre_std_2 =~ s/.fits//;
#($pre_std_2,$post)=split(/\./,$std_2);
    
    $call="radial_sum_rss.pl ".$pre_std_2.".sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 std_2.rss.fits";
    mycall($call);
    $call="img2spec.pl std_2.rss.fits 1 spec_std_2.txt";
    mycall($call);
    
    open(EXPTIME,"read_img_header.pl $std_2 EXPTIME,AIRMASS,DATE-OBS |");
    $line=<EXPTIME>;
    ($file,$exptime,$airmass,$date_obs)=split(" ",$line);
    close(EXPTIME);
#print "$name_std_2 $exptime $airmass\n";
    $factor=1/($exptime); # Filling-factor!
    $date_obs =~ s/T//;
    $date_obs =~ s/\-//g;
    $date_obs =~ s/\://g;
    $extinction=0.15;
    $NC=0;
	while (($date[$NC]<$date_obs)&&($NC<$nc)) {
#    while ($date[$NC]<$date_obs) {
	if ($av[$NC]>$extinction) {
	    $extinction=$av[$NC];
	}
    $NC++;
    }


    $call="calib_det.pl spec_std_2.txt ".$DIR_CONF."/f".$name_std_2.".dat ratio_std2.txt ".$factor." ".$extinction." ".$airmass." 31 15 0 3800 5000";
    mycall($call);
	
    }

}
}

if ($CALIB==1) {
    $call="cp ".$DIR_CONF."/ratio_V1200.txt ratio.txt";
    mycall($call);
}




#




#
# Now we reduce the objects
#

# We count the number of complete object observations:
$call="cat object.txt | grep ' obj' | grep 'p3' > complete.txt";
#$call="cat object.txt | grep ' obj'  > complete.txt";
mycall($call);

open(FH,"wc -l complete.txt |");
$line=<FH>;
($nOBJ,$file)=split(" ",$line);
close(FH);
$nobj=int($nOBJ/3);
print "$nobj completed\n";

$k=0;
open(FH,"<complete.txt");
while($line=<FH>) {
    $line=<FH>;
    $line=<FH>;
    chop($line);
    @data=split(" ",$line);
    $name_obj[$k]=$data[1];
    @parts=split(/\_/,$data[1]);
    $cut="obj_".$parts[1]."_p3_";
    $name_obj[$k] =~ s/$cut//;
    $name_obj[$k] =~ s/ //g;
    $J_obj[$k]=$parts[1];
    #
    # We check that there are the other two pointings!
    # 
    $NAME_JUST=$name_obj[$k];
    open(P1,"cat object.txt | grep ' obj' | grep '$NAME_JUST' | grep 'p1' | wc -l |");
    $NP1=<P1>;
    chop($NP1);
    close(P1);
    open(P2,"cat object.txt | grep ' obj' | grep '$NAME_JUST' | grep 'p2' | wc -l |");
    $NP2=<P2>;
    chop($NP2);
    close(P2);
#    print "$NAME_JUST $NP1 $NP2\n";
    if (($NP1>0)&&($NP2>0)) {
	$k++;
    }
}
close(FH);

$nobj=$k;

#print "nobj=$nobj, k=$k\n";
#exit;



for ($j=0;$j<$nobj;$j++) {
    #
    # We combine the files!
    # 
#    $J=$j+1;
    $J=$J_obj[$j];
    print "Reducing $name_obj[$j] data\n\n";
    $call="cat object.txt | grep ' obj_".$J."_p1' > obj_p1.txt";
    mycall($call);
    $call="cat object.txt | grep ' obj_".$J."_p2' > obj_p2.txt";
    mycall($call);
    $call="cat object.txt | grep ' obj_".$J."_p3' > obj_p3.txt";
    mycall($call);
    $call="cat object.txt | grep ' arc_obj_".$J."' > arc_obj.txt";
    mycall($call);

#
# Identify ARC
#
    open(ARC_OBJ,"<arc_obj.txt");
    while($line=<ARC_OBJ>) {
	chop($line);
	@data=split(" ",$line);
	$arc_obj=$data[0];
    }
    close(ARC_OBJ);

    $call="cat object.txt | grep ' cont_obj_".$J."' > cont_obj.txt";
    mycall($call);
    open(CONT_OBJ,"<cont_obj.txt");
    while($line=<CONT_OBJ>) {
	chop($line);
	@data=split(" ",$line);
	$cont_obj=$data[0];
    }
    close(CONT_OBJ);

    $exptime=1;

    $call="reduce_all.pl ".$arc_obj." ".$cont_obj." ".$arc_obj." ".$exptime." ".$DIR_CONF."/reduce.arc.V1200.config ".$airmass." ".$extinction." > reduce_junk.log";
    mycall($call);

#    $call="split_ppak.pl arc.disp_cor.fits  arc.obj.fits arc.cal.fits arc.sky.fits 0";
#    mycall($call);
    

    $call="Mosaic_rss.pl ".$DIR_CONF."/mos_new.arc.config mos.rss.fits mos.pt.txt";
    mycall($call);
    $call="rss2rss_mask_ppak.pl mos.rss.fits mos.pt.txt mask.rss.fits ".$name_obj[$j]."_arc.rss.fits 1";
    mycall($call);
    $call="rm -r ".$name_obj[$j]."_arc.cube.fits ";
    mycall($call);
    $call="rss2cube_int.pl ".$name_obj[$j]."_arc.rss.fits mos.pt.txt 1 ".$name_obj[$j]."_arc.cube.fits 2 5 1";
    mycall($call);
    $call="cp mos.rss.fits mos_".$name_obj[$j]."_arc.rss.fits";
    mycall($call);
    $call="cp mos.pt.txt mos_".$name_obj[$j]."_arc.pt.txt";
    mycall($call);



}
exit;



sub mycall {
    my $call=$_[0];
    print "$call | ";
    system($call);
    print LOG "$call\n";
    print "DONE\n";
}

sub dumplog {
    my $name=$_[0];
    my $file=$_[1];
    open(DUMP,"<$file");
    my $dd;
    while($dd=<DUMP>) {
	print LOG "# $name | $dd\n";
    }
    close(DUMP);
}


