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
    print "USE: auto_redu.pl DATE(YYYYMMDD)\n";
    exit;
}
$TODAY=$ARGV[0];
$TOMORROW=$TODAY+1;

$DIR_CONF="/disk-b/sanchez/ppak/legacy/";

open(LOG,">auto_redu.log");
open(PROD,">auto_redu.out");


$call="rm -f *cube.fits";
mycall($call);
$call="rm -f *sobj.fits";
mycall($call);
$call="rm -f *rss.fits";
mycall($call);
$call="rm -f *.V.fits";
mycall($call);

$call="cp ".$DIR_CONF."/CCDFlat.fits .";
mycall($call);
$call="cp ".$DIR_CONF."/ARC.dist.id .";
mycall($call);
$call="cp ".$DIR_CONF."/ARC.disp.id .";
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


$call="get_object.pl > object.txt";
mycall($call);


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



if ($nbias>3) {
    $call="imcombine.pl bias.out bias.fits";
    mycall($call);
} else {
    if ($nbias>0) {
	$call="cp ".$file." bias.fits";
	mycall($call);
    }
}
if ($nbias>0) {
    $call="med2df.pl bias.fits mbias.fits 5 5";
    mycall($call);
} else {
#
# If there is no bias, we use a master-one
#
    $call="cp ".$DIR_CONF."/mbias.fits .";
    mycall($call);
}
#
# We combine the SKYFLATS
#
$call="cat object.txt | grep ' sky_1' > sky_1.txt";
mycall($call);
$call="cat object.txt | grep ' sky_2' > sky_2.txt";
mycall($call);

open(SKY_1,"<sky_1.txt");
open(OUT,">sky_1.out");
while($line=<SKY_1>) {
    chop($line);
    @data=split(" ",$line);
    $file=$data[0];
    open(IMSTAT,"imstats.pl $file 350 1014 50 1 | ");
    while($imstat=<IMSTAT>) {
	#print "$file $imstat";
	chop($imstat);
	@dat_imstat=split(" ",$imstat);
	$mean=$dat_imstat[0];
	$max=$dat_imstat[4];
	if (($mean>5000)&&($max<40000)) { 
	    $call="rm -f b".$file;
	    mycall($call);
	    $call="imarith ".$file." - mbias.fits b".$file;
	    mycall($call);
	    print OUT "b$file\n";
	}
    }
    close(IMSTAT);
}
close(OUT);
close(SKY_1);

$call="wc -l sky_1.out > nsky_1.txt";
mycall($call);

open(SKY_1,"<nsky_1.txt");
$line=<SKY_1>;
chop($line);
@data=split(" ",$line);
$n_sky_1=$data[0];
close(SKY_1);



if ($n_sky_1>3) {
    $call="imcombine.pl sky_1.out bsky_1.fits";
    mycall($call);
    $call="rm -f sky_1.fits";
    mycall($call);
    $call="imarith bsky_1.fits + mbias.fits sky_1.fits";
    mycall($call);
} else {
    if ($n_sky_1>0) {
	$call="cp ".$file." sky_1.fits";
	mycall($call);
    }
}


if ($n_sky_1>0) {

#
# Now we reduce the skyflats to produce a fiber-flat
#
    
    $call="cat object.txt | grep cont_sky_1 > cont_sky_1.txt";
    mycall($call);
    open(FH,"<cont_sky_1.txt");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$cont_sky_1=$data[0];
    }
    close(FH);
    
    
    $call="cat object.txt | grep arc_sky_1 > arc_sky_1.txt";
    mycall($call);
    open(FH,"<arc_sky_1.txt");
    while($line=<FH>) {
	chop($line);
	@data=split(" ",$line);
	$arc_sky_1=$data[0];
    }
    close(FH);
    
    $call="reduce_all.pl sky_1.fits ".$cont_sky_1." ".$arc_sky_1."  1 ".$DIR_CONF."/reduce.noflat.config > reduce_junk.log";
    mycall($call);
    $call="fiber_flat_ppak.pl sky_1.disp_cor.fits pfiberflat.fits";
    mycall($call);
    $call="med2df.pl pfiberflat.fits fiberflat.fits 41 1";
    mycall($call);

    
} else {
#
# If there is no skyflat, we use a master-one
#
    $call="cp ".$DIR_CONF."/fiberflat.fits .";
    mycall($call);
    
}

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
    
    $call="reduce_all.pl ".$std_1." ".$cont_std_1." ".$arc_std_1."  1 ".$DIR_CONF."/reduce.flat.config > reduce_junk.log";
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


    $call="calib_det.pl spec_std_1.txt ".$DIR_CONF."/f".$name_std_1.".dat ratio.txt ".$factor." ".$extinction." ".$airmass." 1 11 0 3640 6800";
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
	

	open(EXPTIME,"read_img_header.pl $std_2 EXPTIME,AIRMASS,DATE-OBS |");
	$line=<EXPTIME>;
	($file,$exptime,$airmass,$date_obs)=split(" ",$line);
	close(EXPTIME);
#print "$name_std_2 $exptime $airmass\n";    
	$factor=1/$exptime;
	$date_obs =~ s/T//;
	$date_obs =~ s/\-//g;
	$date_obs =~ s/\://g;
	$extinction=0.15;
	$NC=0;
#	print "PASO $std_2\n";


	while (($date[$NC]<$date_obs)&&($NC<$nc)) {
	    if ($av[$NC]>$extinction) {
		$extinction=$av[$NC];
	    }
#	    print "$NC '$date[$NC]'<$date_obs $av[$NC] $extinction\n";
#	    <stdin>;
	    $NC++;
	}
	
#	print "PASO $std_2\n";
	
	
	
	$call="reduce_all.pl ".$std_2." ".$cont_std_2." ".$arc_std_2."  ".$exptime." ".$DIR_CONF."/reduce.config  ".$airmass." ".$extinction." > reduce_junk.log";
	mycall($call);
    
	$pre_std_2=$std_2;
	$pre_std_2 =~ s/.fits//;
#($pre_std_2,$post)=split(/\./,$std_2);
	
	$call="radial_sum_rss.pl ".$pre_std_2.".sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 std_2.rss.fits";
	mycall($call);
	$call="img2spec.pl std_2.rss.fits 1 spec_std_2.txt";
	mycall($call);
	
	open(CAL,">cal.list");
	$cal_file=$DIR_CONF."/f".$name_std_2.".dat";
	print CAL "$cal_file\n";
	print CAL "$spec_std_2.txt\n";
	close(CAL);
    
	$call="flux_ratio.pl spec_std_2.txt ".$DIR_CONF."/f".$name_std_2.".dat junk_ratio.txt 1 1 0 1 0 10000";
	mycall($call);
	$call="spec_plot.pl junk_ratio.txt calib_test.ps/CPS 0.4 1.6 3650 6800 > calib_test.out";
	mycall($call);
	
	open(FH,"<calib_test.out");
	$line=<FH>;
	close(FH);
	print LOG "# CALIBRATIION ACCURACY $line";
    }

} else {
#
# If there is no calibration star, we use a master-one
#
    $call="cp ".$DIR_CONF."/ratio.txt .";
    mycall($call);
}
#exit;

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
    $call="imcombine.pl obj_p1.txt obj_p1.fits";
    mycall($call);
    $call="imcombine.pl obj_p2.txt obj_p2.fits";
    mycall($call);
    $call="imcombine.pl obj_p3.txt obj_p3.fits";
    mycall($call);

    $call="cat object.txt | grep ' cont_obj_".$J."' > cont_obj.txt";
    mycall($call);
    open(CONT_OBJ,"<cont_obj.txt");
    while($line=<CONT_OBJ>) {
	chop($line);
	@data=split(" ",$line);
	$cont_obj=$data[0];
    }
    close(CONT_OBJ);

    $call="cat object.txt | grep ' arc_obj_".$J."' > arc_obj.txt";
    mycall($call);
    open(ARC_OBJ,"<arc_obj.txt");
    while($line=<ARC_OBJ>) {
	chop($line);
	@data=split(" ",$line);
	$arc_obj=$data[0];
    }
    close(ARC_OBJ);


    open(EXPTIME,"read_img_header.pl obj_p1.fits EXPTIME,AIRMASS,DATE-OBS |");
    $line=<EXPTIME>;
    ($file,$exptime,$airmass,$date_obs)=split(" ",$line);
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

    close(EXPTIME);    
    $exptime=$exptime*5.64;

    $call="reduce_all.pl obj_p1.fits ".$cont_obj." ".$arc_obj." ".$exptime." ".$DIR_CONF."/reduce.config ".$airmass." ".$extinction." > reduce_junk.log";
    mycall($call);

    open(EXPTIME,"read_img_header.pl obj_p2.fits EXPTIME,AIRMASS |");
    $line=<EXPTIME>;
    ($file,$exptime,$airmass,$date_obs)=split(" ",$line);
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

    close(EXPTIME);    
    $exptime=$exptime*5.64;
    $call="reduce_all.pl obj_p2.fits ".$cont_obj." ".$arc_obj." ".$exptime." ".$DIR_CONF."/reduce.config ".$airmass." ".$extinction." > reduce_junk.log";
    mycall($call);

    open(EXPTIME,"read_img_header.pl obj_p3.fits EXPTIME,AIRMASS |");
    $line=<EXPTIME>;
    ($file,$exptime,$airmass,$date_obs)=split(" ",$line);
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

    close(EXPTIME);    
#
# We correct for the fibers aperture
#
    $exptime=$exptime*5.64;

    $call="reduce_all.pl obj_p3.fits ".$cont_obj." ".$arc_obj." ".$exptime." ".$DIR_CONF."/reduce.config ".$airmass." ".$extinction." > reduce_junk.log";
    mycall($call);

#
# If there is a slight change in the transparence, the pointings need to be rescaled
# 

    $call="radial_sum_rss.pl obj_p1.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_1.rss.fits";
    mycall($call);
    $call="img2spec.pl rad_1.rss.fits 4 rad_1.txt";
    mycall($call);
    open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_1.txt 0 |");
    while($L=<FH>) {
	$flux1=$L;
	chop($flux1);
    }
    close(FH);


    $call="radial_sum_rss.pl obj_p2.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_2.rss.fits";
    mycall($call);
    $call="img2spec.pl rad_2.rss.fits 4 rad_2.txt";
    mycall($call);
    open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_2.txt 0 |");
    while($L=<FH>) {
	$flux2=$L;
	chop($flux2);
    }
    close(FH);


    $call="radial_sum_rss.pl obj_p3.sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_3.rss.fits";
    mycall($call);
    $call="img2spec.pl rad_3.rss.fits 4 rad_3.txt";
    mycall($call);
    open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_3.txt 0 |");
    while($L=<FH>) {
	$flux3=$L;
	chop($flux3);
    }
    close(FH);

    print "FLUXES = $flux1 $flux2 $flux3\n";

    if (($flux1>0)&&($flux2>0)) {
	$rat12=$flux2/$flux1;
	if (($rat12<0.95)||($rat12>1.05)) {
	    $call="imarith.pl obj_p2.sobj.fits / ".$rat12." obj_p2.sobj.fits";
	    mycall($call);
	}
    }

    if (($flux1>0)&&($flux3>0)) {
	$rat13=$flux3/$flux1;
	if (($rat13<0.95)||($rat13>1.05)) {
	    $call="imarith.pl obj_p3.sobj.fits / ".$rat13." obj_p3.sobj.fits";
	    mycall($call);
	}
    }


#    exit;

    $call="Mosaic_rss.pl ".$DIR_CONF."/mos.config mos.rss.fits mos.pt.txt";
    mycall($call);
    $call="rm -r ".$name_obj[$j].".cube.fits ";
    mycall($call);
    $call="rss2cube.tcl mos.rss.fits mos.pt.txt ".$name_obj[$j].".cube.fits 3 1e-12 1";
#    $call="rss2cube.tcl mos.rss.fits mos.pt.txt ".$name_obj[$j].".cube.fits 2 1e-12 1";
    mycall($call);
    print PROD "$name_obj[$j].cube.fits\n";

#    exit;

#
# We look for the peak intensity!
#
    $call="get_slice.pl ".$name_obj[$j].".cube.fits img /disk-b/sanchez/ppak/legacy/slice_V.conf";
    mycall($call);
    $img=rfits("img_V_4500_5500.fits");
    ($nx,$ny)=$img->dims;
    $val_max=-1e12;
    $XC=37;
    $YC=37;
    for ($ii=30;$ii<45;$ii++) {
	for ($jj=30;$jj<45;$jj++) {
	    $val=$img->at($ii,$jj);
	    if ($val>$val_max) {
		$val_max=$val;
		$XC=$ii;
		$YC=$jj;
	    }
	}
    }

    $call="cp img_V_4500_5500.fits ".$name_obj[$j].".V.fits";
    mycall($call);

    $call="rm -r ".$name_obj[$j].".scube.fits ";
    mycall($call);
    $call="DAR_det_cube ".$name_obj[$j].".cube.fits  ".$XC." ".$YC." 4 4 ".$name_obj[$j].".scube.fits 20 1050 3 1 0  0 1";
    mycall($call);



    print PROD "$name_obj[$j].scube.fits\n";
    
}

#
# We check if there are partial observations!
#
$call="cat object.txt | grep ' obj'  > total.txt";
mycall($call);
open(FH,"wc -l total.txt |");
$line=<FH>;
($nTOT,$file)=split(" ",$line);
close(FH);
$ntot=int($nTOT/9);
#$nmore=($nTOT-$ntot*9)/3;


for ($i=0;$i<$n_unc;$i++) {
    
#    $J=$nobj+1;
    $J=$Jobj[$i];

    $N_POINTS=0;
    open(FH,"<total.txt");
    while($line=<FH>) {
	chop($line);
	$OBJ_J="obj_".$J;
	if ($line =~ $OBJ_J) {
	    @data=split(" ",$line);
	    $name_obj_now=$data[1];
	    $cut="obj_".$J."_p1_";
	    $name_obj_now =~ s/$cut//;
	    $cut="obj_".$J."_p2_";
	    $name_obj_now =~ s/$cut//;
	    $cut="obj_".$J."_p3_";
	    $name_obj_now =~ s/$cut//;
	    $name_obj_now =~ s/ //g;
	    $N_POINTS++;
	}
    }
    close(FH);
    $nmore=int($N_POINTS/3);
    if ($N_POINTS != 3*$nmore) {
	$nmore++;
    }

    $call="cat object.txt | grep ' cont_obj_".$J."' > cont_obj.txt";
    mycall($call);
    open(CONT_OBJ,"<cont_obj.txt");
    while($line=<CONT_OBJ>) {
	chop($line);
	@data=split(" ",$line);
	$cont_obj=$data[0];
    }
    close(CONT_OBJ);

    $call="cat object.txt | grep ' arc_obj_".$J."' > arc_obj.txt";
    mycall($call);
    open(ARC_OBJ,"<arc_obj.txt");
    while($line=<ARC_OBJ>) {
	chop($line);
	@data=split(" ",$line);
	$arc_obj=$data[0];
    }
    close(ARC_OBJ);    


    print "Reducing $name_obj_now data $nmore\n\n";
#    exit;
    for ($j=0;$j<$nmore;$j++) {
	$n_p_now=0;
	for ($JJ=1;$JJ<4;$JJ++) {
	    $check="p".$JJ."_".$name_obj_now;
	    open(POINT,"cat total.txt | grep '$check' | wc -l |");
	    $point=<POINT>;
	    chop($point);
	    close(POINT);
	    if ($point>0) {
		$p_now[$n_p_now]=$JJ;
		$n_p_now++;
	    }
	}

#	$K=$j+1;
	$K=$p_now[$j];
	$call="cat object.txt | grep ' obj_".$J."_p".$K."' > obj_p".$K.".txt";
	mycall($call);
	$call="imcombine.pl obj_p".$K.".txt obj_p".$K.".fits";
	mycall($call);

	open(EXPTIME,"read_img_header.pl obj_p$K.fits EXPTIME,AIRMASS |");
	$line=<EXPTIME>;
	($file,$exptime,$airmass,$date_obs)=split(" ",$line);
	$date_obs =~ s/T//;
	$date_obs =~ s/\-//g;
	$date_obs =~ s/\://g;

	$extinction=0.15;
	$NC=0;
	while (($date[$NC]<$date_obs)&&($NC<$nc)) {
#	while ($date[$NC]<$date_obs) {
	    if ($av[$NC]>$extinction) {
		$extinction=$av[$NC];
	    }
	    $NC++;
	}

	close(EXPTIME);    
	$call="reduce_all.pl obj_p".$K.".fits ".$cont_obj." ".$arc_obj." ".$exptime." ".$DIR_CONF."/reduce.config ".$airmass." ".$extinction." > reduce_junk.log";
	mycall($call);
	$call="cp obj_p".$K.".sobj.fits ".$name_obj_now."_p".$K.".sobj.fits";
	mycall($call);
	$name_out=$name_obj_now."_p".$K.".sobj.fits";
	#print PROD "$name_obj_now_p$K.sobj.fits\n";
	print PROD "$name_out\n";
    }

    if ($nmore==1) {
	$call="rm -r ".$name_obj_now.".cube.fits ";
	mycall($call);
     
	$call="rss2cube.tcl ".$name_obj_now."_p".$p_now[0].".sobj.fits /disk-b/sanchez/ppak/legacy/ppak_pt_arc.txt ".$name_obj_now.".1P.cube.fits 3 1e-12 1";
	mycall($call);
	print PROD "$name_obj_now.1P.cube.fits\n";
	$call="get_slice.pl ".$name_obj_now.".1P.cube.fits img /disk-b/sanchez/ppak/legacy/slice_V.conf";
	mycall($call);
	$call="cp img_V_4500_5500.fits ".$name_obj_now.".1P.V.fits";
	mycall($call);

	
    } else {


#	exit;


	$call="cp ".$name_obj_now."_p".$p_now[0].".sobj.fits obj_p".$p_now[0].".sobj.fits";
	mycall($call);
	$call="cp ".$name_obj_now."_p".$p_now[1].".sobj.fits obj_p".$p_now[1].".sobj.fits";
	mycall($call);


	$call="radial_sum_rss.pl obj_p".$p_now[0].".sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_1.rss.fits";
	mycall($call);
	$call="img2spec.pl rad_1.rss.fits 4 rad_1.txt";
	mycall($call);
	open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_1.txt 0 |");

    while($L=<FH>) {
	$flux1=$L;
	chop($flux1);
    }
    close(FH);

	
	$call="radial_sum_rss.pl obj_p".$p_now[1].".sobj.fits ".$DIR_CONF."/ppak_pt_arc.txt 5 0 0 rad_2.rss.fits";
	mycall($call);
	$call="img2spec.pl rad_2.rss.fits 4 rad_2.txt";
	mycall($call);
	open(FH,"flux_filter.pl /disk-b/sanchez/ppak/legacy/V_Johnson.txt rad_2.txt 0 |");

    while($L=<FH>) {
	$flux2=$L;
	chop($flux2);
    }
    close(FH);


	print "FLUXES = $flux1 $flux2 \n";
    
	$rat12=$flux2/$flux1;
	if (($rat12<0.95)||($rat12>1.05)) {
	    $call="imarith.pl obj_p".$p_now[1].".sobj.fits / ".$rat12." obj_p".$p_now[1].".sobj.fits";
	    mycall($call);
	}

	$PPP=$p_now[0].$p_now[1];
	$call="Mosaic_rss.pl ".$DIR_CONF."/mos_".$PPP."p.config ".$name_obj_now."_2P.rss.fits mos_".$PPP."p.pt.txt";
	mycall($call);
	$call=$name_obj_now."_2P.rss.fits";
	print PROD "$call\n";
	$call="cp mos_".$PPP."p.pt.txt ".$DIR_CONF."/";
	mycall($call);
	$call="rm -r ".$name_obj_now.".cube.fits ";
	mycall($call);
	$call="rss2cube.tcl ".$name_obj_now."_2P.rss.fits mos_".$PPP."p.pt.txt ".$name_obj_now.".2P.cube.fits 3 1e-12 1";
	mycall($call);
	print PROD "$name_obj_now.2P.cube.fits\n";
	$call="get_slice.pl ".$name_obj_now.".2P.cube.fits img /disk-b/sanchez/ppak/legacy/slice_V.conf";
	mycall($call);
	$call="cp img_V_4500_5500.fits ".$name_obj_now.".2P.V.fits";
	mycall($call);
    }

}

close(PROD);
close(LOG);

#
# We clean the directory
#
$call="rm -f junk*.fits";
mycall($call);
$call="rm -f arc*.fits";
mycall($call);
$call="rm -f sky*.fits";
mycall($call);
$call="rm -f obj_p1.sobj.fits";
mycall($call);
$call="rm -f obj_p2.sobj.fits";
mycall($call);
$call="rm -f obj_p3.sobj.fits";
mycall($call);
$call="rm -f run*b.sobj.fits";
mycall($call);



exit;


sub mycall {
    my $call=$_[0];
    print "$call | ";
    system($call);
    print LOG "$call\n";
    print "DONE\n";
}

