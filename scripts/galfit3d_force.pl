#!/usr/bin/perl
#
#
# This program find peaks in a 2D fiber based spectral image
#
#

#use Statistics::OLS;
#use Math::FFT;
#use Math::Stat;
#use Math::Spline qw(spline linsearch binsearch);
#use Math::Derivative qw(Derivative2);
#use Math::Approx;
#use Astro::FITS::CFITSIO qw( :longnames :constants );
use PDL;
use PDL::Fit::Polynomial; 
use PDL::Filter::Linear;
use PGPLOT;  # Load PGPLOT module
use PDL::Fit::Gaussian;


$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); # "/home/sanchez/sda1/perl/MY/my.pl";

$galfit="/home/sanchez/sda1/galfit/galfit";
$cl="/home/sanchez/sda1/iraf/iraf/unix/hlib/cl.csh";
@param=("name","x","y","mag","re","nsersic","ab","pa");
if ($#ARGV<7) {
    print "USE: galfit3d_force.pl INPUT_CUBE.fits MOD_CUBE.fits RES_CUBE.fits GALFIT_CONF INPUT_LOGFILE INPUT_NPOLY OUTPUT_LOGFILE QUALITY.txt [PLOT]\n";
    print "INPUT_NPOLY: File in the form\n";
    print "Model_number parameter_number Poly_Order_to_FIT_and_FIX\n"; 
    print "Ify Poly_Order_to_FIT_and_FIX it is not an integer, \n";
    print "then sustitute for this value.\n";
    print "Models numbering starts with 0\n";
    print "Parameters are: @param\n";
    print "GALFIT_CONF: A) object.fits, B) output.fits\n";
    exit;
}

$input_cube=$ARGV[0];
$mod_cube=$ARGV[1];
$res_cube=$ARGV[2];
$galfit_conf=$ARGV[3];
$in_logfile=$ARGV[4];
$in_polyfile=$ARGV[5];
$out_logfile=$ARGV[6];
$qfile=$ARGV[7];
$plot=0;
if ($#ARGV==8) {
    $plot=1;
}

#
# We create a tmp galfit:
#
print "Reading $galfit_conf\n";
$model_number=0;
$nh=0;
open(FH,"<$galfit_conf");
while($line=<FH>) {
    chop($line);
#    if ($line !~ "#") {
	if ($nh==0) {
	    do {	
#	    print "$line\n $nh\n";
		$header[$nh]=$line;
		$line=<FH>;
		chop($line);
		$nh++;
	    } while ($line !~ /S\)/);
	    $header[$nh]=$line;
	    $line=<FH>;
	    chop($line);
	    $nh++;
	}
	
	@data=split(" ",$line);
#	print "$line\ndata[0]=$data[0]\n";
	if ($data[0] =~ /1\)/) {
#	@data=split(" ",$line);
	    $a_fit[$model_number][1]=$data[3];
	    $a_fit[$model_number][2]=$data[4];
	}
	if ($data[0] =~ /3\)/) {
#	@data=split(" ",$line);
	    $a_fit[$model_number][3]=$data[2];
	}
	if ($data[0] =~ /4\)/) {
#	@data=split(" ",$line);
	    $a_fit[$model_number][4]=$data[2];
#	    print "$line\n$model_number $a_fit[$model_number][4]=$data[2];\n";
	}
	if ($data[0] =~ /5\)/) {
#	@data=split(" ",$line);
	    $a_fit[$model_number][5]=$data[2];
	}
	if ($data[0] =~ /8\)/) {
#	@data=split(" ",$line);
	    $a_fit[$model_number][6]=$data[2];
	}
	if ($data[0] =~ /9\)/) {
#	@data=split(" ",$line);
	    $a_fit[$model_number][7]=$data[2];
	}
#	print "$data[2]\n";
	if ($data[0] =~ /Z\)/) {
#	    @data=split(" ",$line);
	    $a_subs[$model_number]=$data[1];
	    for ($j=0;$j<8;$j++) {
		$a_fix[$model_number][$j]=-666;
	    }
#
#
#
#	    for ($j=0;$j<8;$j++) {
#		print "a_fit[$model_number][$j]=$a_fit[$model_number][$j]\n"
#		}
	    
	    
	    $model_number++;		
	}
#}	
}
close(FH);
print "Number of Models=$model_number\n";
#exit;
#
# We check for the fixing paramters:
#
$n_fix=0;
open(FIT,"<$in_polyfile");
while($line=<FIT>) {
    chop($line);
    if ($line !~ "#") {
	@array=split(" ",$line);
	$n_mod=$array[0];
	$n_par=$array[1];
	$a_fix[$n_mod][$n_par]=$array[2];
	$a_box[$n_mod][$n_par]=$array[3];
	$n_fix++;
    }
}
close(FIT);

#exit;

#
# We read the input log_file
#


print "Reading input GALFIT log file: $in_log_file\n";
open(LOG,"<$in_logfile");
while ($line=<LOG>) {
    if ($line =~ "Output") {
	$line=<LOG>; #blank line
	for ($i=0;$i<$model_number;$i++) {
	    $results=<LOG>;
	    $errors=<LOG>;
	    $results =~ s/\://;
	    $errors =~ s/\://;
	    $results =~ s/\,//;
	    $errors =~ s/\,//;
	    $results =~ s/\(//;
	    $errors =~ s/\(//;
	    $results =~ s/\)//;
	    $errors =~ s/\)//;
	    @data=split(" ",$results);
	    for ($j=0;$j<8;$j++) {
		$a_param[$n][$i][$j]=$data[$j];
	    }
	    @data=split(" ",$errors);
	    for ($j=0;$j<8;$j++) {
		$a_e_param[$n][$i][$j]=$data[$j];
	    }	
	}
	$n++;
    }
}
close(LOG);
print "$n index readed. DONE\n";

#
# We perform the polynomical fittings!
#
for ($i=0;$i<$model_number;$i++) {
    for ($j=1;$j<8;$j++) {    
	if ($a_fix[$i][$j]!=-666) {
	    if ($a_fix[$i][$j]==int($a_fix[$i][$j])) {		
		print "We fit parameter '$param[$j]' of '$a_param[0][$i][0]' with a polynomical function of order '$a_fix[$i][$j]'\n";
		for ($k=0;$k<$n;$k++) {
		    $a_x[$k]=$k;
		    $a_y[$k]=$a_param[$k][$i][$j];
#		    print "$k $a_y[$k]\n";
		}
		$med_box=$a_box[$i][$j];
		if ($med_box>3) {
		   my @tmp=median_filter($med_box,\@a_y); 
		   @a_y=@tmp;
		}

		$npoly=$a_fix[$i][$j];
		my $y_pdl = pdl(@a_y);
		($s_y,$coeff) = fitpoly1d $y_pdl,$npoly;
		@a_y_out=list($s_y);
		($min, $max) = minmax(@a_y);
		if ($plot==1) {
		    pgbegin(0,"/xs",1,1);
		    pgsfs(1.2);
		    pgscf(2);             # Set character font
		    pgslw(2);             # Set line width
		    pgsch(1.6);           # Set character height
		    pgenv(1,$n,$min-0.1,$max+0.1,0,0);
		    pglabel("Index","$param[$j]","n.poly=$a_fix[$i][$j]");
		    pgsci(2);
		    pgpoint($n,\@a_x,\@a_y,1);
		    pgsci(1);
		    pgline($n,\@a_x,\@a_y_out);
		    pgclose;
		    pgend;
		    print "Press Enter\n";
		    <stdin>;
		}
		for ($k=0;$k<$n;$k++) {
		    $a_param[$k][$i][$j]=$a_y_out[$k];
		}				
	    } else {
		print "We fix parameter '$param[$j]' of '$a_param[0][$i][0]' to the value '$a_fix[$i][$j]'\n";
		for ($k=0;$k<$n;$k++) {
		    $a_param[$k][$i][$j]=$a_fix[$i][$j];
		}
	    }
	}
    }
}
#exit;
if ($plot==1) {
    print "Press Enter to start Fitting\n";
    <stdin>;
}

#$pdl=rfits($input_cube,{data=>0});
#$nx=$pdl->{"NAXIS1"};
#$ny=$pdl->{"NAXIS2"};
#$nz=$pdl->{"NAXIS3"};



$pdl_input_cube=rfits($input_cube);
$nx=$pdl_input_cube->hdr->{"NAXIS1"};
$ny=$pdl_input_cube->hdr->{"NAXIS2"};
$nz=$pdl_input_cube->hdr->{"NAXIS3"};
$pdl_mod_cube=rfits($input_cube);
$pdl_res_cube=rfits($input_cube);


system("rm fit.log");
#
#We create the output files
#
#system("rm $mod_cube");
#system("rm $res_cube");
#open(TMP,">tmp.cl");
#print TMP "\n";
#print TMP "\n";
#print TMP "\n";
#print TMP "imcopy $input_cube $mod_cube\n";
#print TMP "imcopy $input_cube $res_cube\n";
#print TMP "\n";
#print TMP "logout\n";
#print TMP "\n";
#close(TMP);
#$call=$cl." < tmp.cl > junk.tmp";
#system($call);
 


#exit;


#
# Loop over the Z direction:
#
for ($k=0;$k<$nz;$k++) {
    $kk=$k+1;
    print "$k/$nz...";
    my $slice=$pdl_input_cube->slice(",,($k)");
#    $slice->wfits("object.fits");

    system("rm object.fits");
    open(TMP,">copy.cl");
    print TMP "\n";
    print TMP "imcopy $input_cube\[*,*,$kk\] object.fits\n";
    print TMP "\n";
    print TMP "logout\n";
    print TMP "\n";
    close(TMP);
    $call=$cl." < copy.cl > junk.tmp";
    system($call);
    print "copied...";

#
# We create a tmp galfit:
#

    open(FH,">galfit.tmp");
    for ($i=0;$i<$nh;$i++) {
	print FH "$header[$i]\n";
    }
    for ($i=0;$i<$model_number;$i++) {
	$ii=$i+1;
	print FH "\n";
	print FH "#Object number: $ii\n";
	print FH " 0) $a_param[0][$i][0] # Object type\n";
	print FH " 1) $a_param[$k][$i][1] $a_param[$k][$i][2] $a_fit[$i][1] $a_fit[$i][2] # ($param[1],$param[2])\n";
	print FH " 3) $a_param[$k][$i][3] $a_fit[$i][3] # ($param[3])\n";
	print FH " 4) $a_param[$k][$i][4] $a_fit[$i][4] # ($param[4])\n";
	print FH " 5) $a_param[$k][$i][5] $a_fit[$i][5] # ($param[5])\n";
	print FH " 6) 0.000      0       #     ----- \n";
	print FH " 7) 0.000      0       #     ----- \n";
	if ($a_param[$k][$i][6]==0) {
	    $a_param[$k][$i][6]=1;
	}
	print FH " 8) $a_param[$k][$i][6] $a_fit[$i][6] # ($param[6])\n";
	print FH " 9) $a_param[$k][$i][7] $a_fit[$i][7] # ($param[7])\n";
	print FH "10) 0.000      0       #  diskiness(-)/boxiness(+)\n";
	print FH " Z) $a_subs[$i]        #  Output option (0 = residual, 1 = Don't subtract) \n";
#$a_param[$k][$i][$j];
    }
    close(FH);

#    print "BREAK"; exit;
#    $call=$galfit." ".$galfit_conf;
    $call=$galfit." galfit.tmp > galfit.junk";
    system($call);
    system ("rm galfit.??");
    system ("rm galfit.junk");
    print "modelled...";
#    open(TMP,">tmp.cl");
#    print TMP "\n";
#    print TMP "imcopy output.fits\[2\] $mod_cube\[*,*,$kk\]\n";
#    print TMP "\n";
#    print TMP "imcopy output.fits\[3\] $res_cube\[*,*,$kk\]\n";
#    print TMP "\n";
#    print TMP "logout\n";
#    print TMP "\n";
#    print TMP "display output.fits\[1\] 1 zs- zr- z1=-0.001 z2=0.02\n";
#    print TMP "display output.fits\[2\] 2 zs- zr- z1=-0.001 z2=0.02\n";
#    print TMP "display output.fits\[3\] 3 zs- zr- z1=-0.001 z2=0.02\n";
#    close(TMP);
#    $call=$cl." < tmp.cl > junk.tmp";
#    system($call);

    my $mod;
    open(DIR,"ls output.fits |");
    $dir=<DIR>;
    close(DIR);
    chop($dir);
#    print "\n'$dir'\n";
    $fit=0.0;
    if ($dir !~ "output.fits" ) {
	$mod=zeroes($nx,$ny);
    } else {
	$mod=rfits("output.fits[2]");
	$fit=1.0;
    }    
    ($nx_c,$ny_c)=$mod->dims;
#    $mod->wfits("test_mod.fits");
    for ($i=0;$i<$nx_c;$i++) {
	for ($j=0;$j<$ny_c;$j++) {	
	    $val=$mod->at($i,$j);
	    set($pdl_mod_cube,$i,$j,$k,$val);
	}
    }
    my $res;
    if ($dir !~ "output.fits" ) {
	$res=$slice;
    } else {
	$res=rfits("output.fits[3]");
    }    

    for ($i=0;$i<$nx_c;$i++) {
	for ($j=0;$j<$ny_c;$j++) {	
	    $val=$res->at($i,$j);
	    set($pdl_res_cube,$i,$j,$k,$val);
	}
    }


    print "DONE\n";
    open(TMP,">display.cl");
    print TMP "\n";
    print TMP "display output.fits\[1\] 1 zs- zr- z1=-0.001 z2=0.02\n";
    print TMP "display output.fits\[2\] 2 zs- zr- z1=-0.001 z2=0.02\n";
    print TMP "display output.fits\[3\] 3 zs- zr- z1=-0.001 z2=0.02\n";
    close(TMP);
#    $call=$cl." < display.cl > junk.tmp";
#    system($call);
#    print "DONE\n";
#    print "Press Enter"; <stdin>;   
    if ($k==0) {
	open(Q,">$qfile");
    } else {
	open(Q,">>$qfile");
    }
    print Q "$k $fit\n";
    close(Q);   
    system("cp fit.log $out_logfile");
}
print "Model done\n";
$pdl_mod_cube->wfits($mod_cube);
$pdl_res_cube->wfits($res_cube);
print "$mod_cube and $res_cube saved\n";


exit;




