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

if ($#ARGV<2) {
    print "USE: Mosaic_rss_overlap.pl CONFIG.txt RSS.FITS POSITION_TABLE.txt CUT [DIST_MIN] [NS1 NS2]\n";
    print "CONFIG: infile in_position_table dx dy\n";
    exit;
}

$config=$ARGV[0];
$outfile=$ARGV[1];
$outpt=$ARGV[2];
$cut=$ARGV[3];
$DIST_MIN=1;
if ($#ARGV==4) {
    $DIST_MIN=$ARGV[4];
}






$nf=0;
$dxmin=1e12;
$dxmax=-1e12;
$dymin=1e12;
$dymax=-1e12;
open(FH,"<$config");
while($line=<FH>) {
#    chop($line);
    @data=split(" ",$line);
    $infile[$nf]=$data[0];
    $inpt[$nf]=$data[1];
    $dx[$nf]=$data[2];
    $dy[$nf]=$data[3];
    if ($dxmin>$dx[$nf]) {
	$dxmin=$dx[$nf];
    }
    if ($dxmax<$dx[$nf]) {
	$dxmax=$dx[$nf];
    }
    if ($dymin>$dy[$nf]) {
	$dymin=$dy[$nf];
    }
    if ($dymax<$dy[$nf]) {
	$dymax=$dy[$nf];
    }
    $nf++;
}
close(FH);
print "$nf files\n";

#
# Reading Position tables:
#
open(OUT,">$outpt");
$n=0;
for ($j=0;$j<$nf;$j++) {
    $pt_file=$inpt[$j];
    print "Adding $pt_file, ";
    open(FH,"<$pt_file");
    if ($j==0) {
	$line=<FH>;
	print OUT "$line";
    } else {
	$line=<FH>;
    }
    $kk=0;
    while($line=<FH>) {
	chop($line);
        @data=split(" ",$line);
	$id=$data[0]+$j*$n;
	if ($j==0) {
	    $n++;
	}
	$x[$kk]=$data[1]+$dx[$j];
	$y[$kk]=$data[2]+$dy[$j];
	$match[$kk]=-1;
	$kk++;
	$type=$data[3];
#	print OUT "$id $x $y $type\n";
    }
    close(FH);
    print "Position Table readed\n";
    $n_match=0;
    if ($j==0) {
	for ($nk=0;$nk<$kk;$nk++) {
	    $xx[$nk]=$x[$nk];
	    $yy[$nk]=$y[$nk];	    
	  #  print "$xx[$nk] $yy[$nk]\n";
	}
	#<stdin>;
	$nk=$kk;
    } else {
	for ($ii=0;$ii<$nk;$ii++) {
	    for ($i=0;$i<$kk;$i++) {
		$dist=sqrt(($x[$i]-$xx[$ii])**2+($y[$i]-$yy[$ii])**2);
		#print "$j $i ($x[$i],$y[$i]) $ii ($xx[$i],$yy[$ii]) $dist\n";
#		if ($dist<2) {
		if ($dist<$DIST_MIN) {
		    $match[$i]=$ii; # There is a match
		    $n_match++;
#		    print "$i $ii $n_match\n";
		}
	    }
#	    <stdin>;
	}
    }
#    print "N.Match=$n_match , ";
#    <stdin>;
    $a_in=rfits($infile[$j]);
#    $a_test=rfits($infile[$j]);
#    $a_in=$a_test->slice("700:850,:");
    ($nx,$ny)=$a_in->dims;
    $Ns1=int(0.25*$nx);
    $Ns2=int(0.75*$nx);

    if ($#ARGV==6) {
	$DIST_MIN=$ARGV[4];
	$Ns1=$ARGV[5];	
	$Ns2=$ARGV[6];
    }


    print "NS=$Ns1,$Ns2\n";
    if ($j==0) {
	$a_new=$a_in;
	$ny_pre=0;
	$ny_new=$ny;
    } else {
	($nx_pre,$ny_pre)=$a_out->dims;
	$h=$a_in->gethdr;
#	$h=$a_test->gethdr;


	if ($nx<$nx_pre) {
	    $a_tmp=zeroes($nx_pre,$ny);
	    for ($I=0;$I<$nx;$I++) {
		my $t = $a_tmp->slice("$I,:");
		$t .= $a_in->slice("$I,:");
	    }
	    $a_in=$a_tmp;
	}
	($nx,$ny)=$a_in->dims;
	$ny_new=$ny_pre+$ny-$n_match;
	print "($nx,$ny) ($nx_pre,$ny_pre) $ny_new\n";

	$a_new=zeroes($nx_pre,$ny_new);
	for ($i=0;$i<$nx;$i++) {
	    $a_x[$i]=$i;
	}
#
# We add the points with single data
#
	$ny_sec=$ny_pre-1;
	$t=$a_new->slice(":,0:$ny_sec");
	print "PASO---\n";
	$t .=$a_out;

	
	if ($n_match>0) {
	    $nnm=0;
	    my $rat_pdl=ones($nx,$n_match);
	    for ($jj=0;$jj<$ny;$jj++) {
		if ($match[$jj]!=-1) {
		    $ii=$match[$jj];
#		    @DIMS=$rat_pdl_now->dims;
		    for ($i=0;$i<$nx;$i++) {
			$A1=$a_in->at($i,$jj);
			$A2=$a_new->at($i,$ii); 
			if ($A2>$cut) {
			    $val_now=$A1/$A2;
			} else {
			    $A2=1;
			    $val_now=1;
			}
#			$val_now=$rat_pdl_now->at($i);
#			print "PASO $jj/$ny $nnm $i/$nx ($DIMS[0],$DIMS[1])\n";
			set($rat_pdl,$i,$nnm,$val_now);
		    }
		    $nnm++;
		}
	    }
#	    $kernel=ones(10);
#	    $kernel=$kernel/10;
#	    $rat_smooth = med2d $rat_pdl,$kernel;



#	    $ratio=average($rat_pdl->xchg(0,1));
#	    $ratio=average($rat_pdl->xchg(0,1));
	    $ratio=medover($rat_pdl->xchg(0,1));
	    ($meanR,$rmsR,$medianR,$minR,$maxR) = stats($ratio->slice("$Ns1:$Ns2"));
	    print "MEAN=$meanR, RMS=$rmsR $Ns1,$Ns2 $n_match\n";
#	    if ($rmsR>0.5*$meanR) {
#		$ratio=ones($nx)*$meanR;
#	    }
	    $ratio_tmp=$ratio;
	    @list_ratio_tmp=list($ratio_tmp);
	    @list_ratio_med=median_filter(50,\@list_ratio_tmp);
	    for ($iii=0;$iii<30;$iii++) {
		$list_ratio_med[$iii]=$list_ratio_med[31];
	    }

	    for ($iii=$nx-30;$iii<$nx;$iii++) {
		$list_ratio_med[$iii]=$list_ratio_med[$nx-31];
	    }

	    $npoly=4;
	    ($ratio,$coeff) = fitpoly1d(pdl(@list_ratio_med),$npoly);

	    
	    @list_ratio=list($ratio);
	    
	    for ($iii=0;$iii<$nx;$iii++) {
		$NX[$iii]=$iii;
	    }

	    pgbegin(0,"/xs",1,1);
	    pgsfs(1.2);
	    pgscf(2);             # Set character font
	    pgslw(2);             # Set line width
	    pgsch(1.2);           # Set character height
	    pgenv(0,$nx,0,4,0,0);
	    pglabel("Wavelength","Flux","");
	    pgsci(1);
	    pgpoint($nx,\@NX,\@list_ratio_tmp,1);
	    pgsci(8);
	    pgpoint($nx,\@NX,\@list_ratio_med,22);
	    pgsci(2);
	    pgline($nx,\@NX,\@list_ratio);
	    pgclose;
	    pgend;



	} else {
	    $ratio=ones($nx);
	}
	$DIMS=$ratio->dims;
#	print "$DIMS[0] $DIMS[1]\n";

	$cont=0;
	$nnm=0;
	for ($jj=0;$jj<$ny;$jj++) {
	    $jjj=$cont+$ny_pre;
	    if ($match[$jj]==-1) {
		$val=$a_in->slice(":,$jj");
		$val=$val/$ratio;
		$t=$a_new->slice(":,$jjj"); 
		$t .=$val;#->slice(":,0:400"));
		$xx[$nk]=$x[$jj];
		$yy[$nk]=$y[$jj];
		$nk++;
		$cont++;
	    } else {
		$ii=$match[$jj];
		$val=$a_in->slice(":,$jj");
		$val=$val/$ratio;
#		$val=$val/$rat_val[$nnm];
		$t=$a_new->slice(":,$ii"); 
		$t .=($t+$val)/2;
		$nnm++;
	    }
	}
	
    }
    if ($j==0) {
	$sum_ratio=1;
    } else {
	$sum_ratio=average($ratio->slice("$Ns1:$Ns2"));
    }

#    print "$ratio\n";
    print "NY_pre=$ny_pre, NY_new=$ny_new NK=$nk RATIO=$sum_ratio\n";
    $a_out=$a_new;
    $a_new=zeroes(1,1);
    print "file $j added\n";
    ($nx_pre,$ny_pre)=$a_out->dims;
#    print "PASO===($nx_pre,$ny_pre)\n";
}

for ($j=0;$j<$nk;$j++) {
    $id=$j+1;
    print OUT "$id $xx[$j] $yy[$j] 1\n";
}

close(OUT);



#$h=$a_in->gethdr;
#$h = {NAXIS=>2, NAXIS1=>$nx, NAXIS=>$ny_new, COMMENT=>"Sample FITS-style header"};

#
#$a_out->wfits($outfile);
#$a_out->rfits($outfile);
$a_out->sethdr($h);
$a_out->hdr->{NAXIS2}=$ny_new;
$a_out->hdr->{NAXIS1}=$nx_pre;
$a_out->wfits($outfile);

print "DONE\n";
exit;



