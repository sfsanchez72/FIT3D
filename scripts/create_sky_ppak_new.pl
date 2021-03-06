#!/usr/bin/perl
use PGPLOT;  # Load PGPLOT module
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




$F3D=$ENV{FIT3D_PATH}; require("$F3D/scripts/my.pl"); #("/home/sanchez/sda2/code/R3D/my.pl");

if ($#ARGV<4) {
    print "USE: create_sky_ppak_new.pl INPUT OBJECT SKY npoly[NOT USED] PLOT\n";
    exit;
}

$inputfile=$ARGV[0];
$objfile=$ARGV[1];
$skyfile=$ARGV[2];
$npoly=$ARGV[3];
$spy=$ARGV[4];

$n=0;
$nc=15;
$ns=36;
#open(FH,"<ppak_positions.txt");
while($line=<DATA>) {
    chop($line);
    @data=split(" ",$line);
    $type[$n]=0;
    if (($data[1]>400)&&($data[1]<500))    {
	$type[$n]=1;
#	print "\n";
    } 

    if (($data[1]>500)&&($data[1]<900))    {
	$type[$n]=2;
#	print "\n";
    }

    if ($data[1]<400)    {
#	print "$data[0] ";
    }


    if ($data[1]<900) {
	$n++;
    }

    
}
#close(FH);

#print "N=$n\n";

print "Reading input file\n";
$naxes=read_naxes($inputfile);
($n1,$n2) = @$naxes;
@input=read_img($inputfile);
print "Done\n";
print "$n1 $n2\n";


$jc=0;
$js=0;
$ji=0;
print "triggering\n";
for ($j=0;$j<$n2;$j++) {
$y_max=-100000;
$y_min=100000;

#    print "type=$type[$j] $jc/15 $js/36\n";
    if ($type[$j]==2) {
	for ($i=0;$i<$n1;$i++) {
	    $cal[$jc][$i]=$input[$j][$i];
	    $y[$i]=$input[$j][$i];
	    $x[$i]=$i;
	    if ($y_max<$y[$i]) {
		$y_max=$y[$i];
	    }
	    if ($y_min>$y[$i]) {
		$y_min=$y[$i];
	    }
	}


	$jc++;

    }

    if ($type[$j]==1) {
	for ($i=0;$i<$n1;$i++) {
	    $sky[$js][$i]=$input[$j][$i];
	    $pos[$js]=$j;
	    $y[$i]=$input[$j][$i];
	    $x[$i]=$i;
	    if ($y_max<$y[$i]) {
		$y_max=$y[$i];
	    }
	    if ($y_min>$y[$i]) {
		$y_min=$y[$i];
	    }
	}




	    $js++;

    }

    if ($type[$j]==0) {
	for ($i=0;$i<$n1;$i++) {
	    $obj[$ji][$i]=$input[$j][$i];
	    $y[$i]=$input[$j][$i];
	    $x[$i]=$i;
	    if ($y_max<$y[$i]) {
		$y_max=$y[$i];
	    }
	    if ($y_min>$y[$i]) {
		$y_min=$y[$i];
	    }
	}



	$ji++;

    }


#    print "$j $jc $js $type[$j]\n";
}
for ($j=0;$j<382;$j++) {
    $new_pos[$j]=$j;
}
my @sky_331;
for ($i=0;$i<$n1;$i++) {
    my @array;
    $y_min=1e12;
    $y_max=-1e12;
    $kk=0;
    for ($j=0;$j<$js;$j++) {
	$val_now=$sky[$j][$i];
#	if ($val_now>0) {
	    $array[$kk]=$sky[$j][$i];
	    $posn[$kk]=$pos[$j];
	    if ($y_min>$array[$kk]) {
		$y_min=$array[$kk];
	    }
	    if ($y_max<$array[$kk]) {
		$y_max=$array[$kk];
	    }
	    $kk++
    }

    $kk2=int(0.5*$kk);
    for ($jj=$kk2;$jj>-1;$jj--) {
	if ($array[$jj]==0) {
	    $array[$jj]=$array[$jj+1];
	}
    }
    for ($jj=$kk2;$jj<$kk;$jj++) {
	if ($array[$jj]==0) {
	    $array[$jj]=$array[$jj-1];
	}
    }
    my $med=median(@array);
    my $sig=sigma(@array);

    for ($jj=0;$jj<$kk;$jj++) {
	$marray[$jj]=$array[$jj];
	if (abs($marray[$jj]-$med)>2*$sig) {
	    $marray[$jj]=$med;
	}
    }
#    @marray=median_clip1d(\@array,36,2,1);




#    @marray=median_filter(10,\@array);
#    print "$#array $#marray\n";

#    print "PASO\n";


#    ($s_f,$coeff) = fitpoly1d(pdl(@pos),pdl(@marray),$npoly);


#    my $pdl_array = interpol(pdl(@new_pos), pdl(@pos), $s_f);
    my $pdl_array = interpol(pdl(@new_pos), pdl(@pos), pdl(@marray));


    my @new_array=list($pdl_array);
    if ($spy==1) {
	pgbegin(0,"/xs",1,1);
	pgsfs(1.2);
	pgscf(2);             # Set character font
	pgslw(2);             # Set line width
	pgsch(1.6);           # Set character height
	pgenv(0,382,$y_min,$y_max,0,0);
	pgsch(1.6);           # Set character height
	pglabel("Spec.Id","Counts","Row $i");
	pgsch(2.2);           # Set character height
	pgsci(1);
	pgpoint($js,\@pos,\@array,2);
	pgsci(3);
	pgpoint($js,\@pos,\@marray,1);
	pgsci(3);
	pgpoint($js,\@pos,\@marray,1);
	pgsci(4);
	pgpoint($js,\@pos,\@array,2);
#
	pgsci(2);
	pgline(382,\@new_pos,\@new_array);    
	pgsci(1);
	pgclose;
	pgend;
	<stdin>;
    }
    $ji=0;
    for ($j=0;$j<382;$j++) {
	if ($type[$j]==0) {
	    $sky_331[$ji][$i]=$new_array[$j];
	    $ji++;
	}
    }
}

#print "Writting the results\n";
#system("rm $calfile");
#write_img($calfile,$n1,15,\@cal);
#print "$calfile saved\n";
system("rm $skyfile");
write_img($skyfile,$n1,331,\@sky_331);
print "$skyfile saved\n";
system("rm $objfile");
write_img($objfile,$n1,331,\@obj);
print "$objfile saved\n";

exit;


__DATA__
      1     501    501       501    
      2     502    502       502    
      3     999    999       999    
      4       1      0.780     0.000
      5       2     -0.878    -0.169
      6     401    401       401    
      7       3      0.000    -0.676
      8       4      0.780    -0.676
      9       5     -0.780    -0.676
     10       6      0.585    -0.338
     11       7     -0.390    -1.013
     12       8      0.390    -1.013
     13       9      0.975    -0.338
     14      10     -0.488    -0.507
     15      11      0.390    -0.676
     16      12      0.000    -1.013
     17     402    402       402    
     18      13     -0.975    -0.338
     19      14      1.267    -0.169
     20      15     -0.390    -0.676
     21      16      0.683    -0.844
     22      17     -0.683    -0.169
     23      18      0.098    -0.844
     24      19     -0.683    -0.844
     25      20      1.170     0.000
     26      21     -1.267    -0.169
     27      22      0.683    -0.507
     28     403    403       403    
     29      23     -0.683    -1.182
     30      24      0.293    -0.844
     31      25      0.878    -0.844
     32      26     -0.683    -0.507
     33     503    503       503    
     34      27      0.878    -0.169
     35      28     -0.293    -0.844
     36      29      0.098    -1.182
     37      30     -1.072    -0.507
     38     404    404       404    
     39      31      1.072    -0.507
     40      32     -0.585    -0.338
     41      33     -0.293    -1.182
     42      34      0.195    -0.676
     43      35      0.488    -1.182
     44      36     -0.195    -0.676
     45      37     -0.878    -0.844
     46      38      0.683    -0.169
     47      39      0.780    -1.013
     48      40     -0.585    -0.676
     49     405    405       405    
     50      41      1.072    -0.169
     51      42     -0.780    -0.338
     52      43      0.585    -0.676
     53      44     -0.195    -1.013
     54      45     -0.975    -0.676
     55      46      0.683    -1.182
     56      47      1.170    -0.338
     57      48     -0.585    -1.013
     58      49      0.488    -0.844
     59      50     -1.072    -0.169
     60     406    406       406    
     61      51     -0.098    -0.844
     62      52      0.293    -1.182
     63      53      0.878    -0.507
     64      54     -0.780    -1.013
     65     504    504       504    
     66      55     -0.878    -0.507
     67      56      0.975     0.000
     68      57     -0.098    -1.182
     69      58      0.585    -1.013
     70     407    407       407    
     71      59     -1.170    -0.338
     72      60     -0.488    -0.844
     73      61      0.488    -0.507
     74      62      1.365     0.000
     75      63      0.195    -1.013
     76      64      0.975    -0.676
     77      65     -0.488    -1.182
     78      66      0.780    -0.338
     79      67     -1.365    -0.338
     80      68      0.390    -1.351
     81     408    408       408    
     82      69      1.365    -0.338
     83      70     -0.975    -1.013
     84      71     -0.488    -1.520
     85      72      1.072    -0.844
     86      73      1.658    -0.169
     87      74     -1.365    -0.676
     88      75      0.000    -1.351
     89      76      0.780    -1.351
     90      77     -0.975    -1.351
     91      78      1.365    -0.676
     92     409    409       409    
     93      79      0.293    -1.520
     94      80     -1.170    -0.676
     95      81      1.560     0.000
     96      82     -0.390    -1.351
     97     505    505       505    
     98      83     -1.755    -0.338
     99      84      1.072    -1.182
    100      85      1.560    -0.338
    101      86     -1.170    -1.013
    102     410    410       410    
    103      87     -0.195    -1.689
    104      88      0.683    -1.520
    105      89      1.267    -0.844
    106      90     -0.780    -1.351
    107      91     -1.462    -0.507
    108      92      1.950     0.000
    109      93     -1.365    -1.013
    110      94      1.072    -1.520
    111      95     -0.195    -1.351
    112      96     -1.658    -0.169
    113     411    411       411    
    114      97     -0.975    -1.689
    115      98      1.267    -0.507
    116      99      0.098    -1.520
    117     100      1.365    -1.013
    118     101     -1.560    -0.676
    119     102      0.878    -1.182
    120     103     -1.072    -1.182
    121     104      0.780    -1.689
    122     105     -0.683    -1.520
    123     106     -1.462    -0.169
    124     412    412       412    
    125     107      1.755     0.000
    126     108      0.390    -1.689
    127     109      1.658    -0.507
    128     110     -1.267    -0.844
    129     506    506       506    
    130     111     -1.170    -1.351
    131     112      0.975    -1.013
    132     113      0.195    -1.351
    133     114     -0.585    -1.689
    134     413    413       413    
    135     115      1.462    -0.169
    136     116     -1.852    -0.169
    137     117     -0.098    -1.520
    138     118      0.975    -1.351
    139     119     -1.267    -0.507
    140     120      1.462    -0.844
    141     121      0.488    -1.520
    142     122     -1.072    -0.844
    143     123     -0.878    -1.520
    144     124      1.170    -0.676
    145     414    414       414    
    146     125      1.755    -0.338
    147     126      1.170    -1.013
    148     127     -1.560    -0.338
    149     128      0.878    -1.520
    150     129     -1.267    -1.182
    151     130     -0.293    -1.520
    152     131      1.462    -0.507
    153     132      1.170    -1.351
    154     133      0.195    -1.689
    155     134     -1.658    -0.507
    156     415    415       415    
    157     135      1.852    -0.169
    158     136      0.585    -1.351
    159     137     -1.072    -1.520
    160     138      1.560    -0.676
    161     507    507       507    
    162     139     -0.390    -1.689
    163     140      1.267    -1.182
    164     141      0.975    -1.689
    165     142     -0.585    -1.351
    166     416    416       416    
    167     143     -1.462    -0.844
    168     144      0.000    -1.689
    169     145     -0.878    -1.182
    170     146     -0.780    -1.689
    171     147      0.585    -1.689
    172     148     -0.585     0.000
    173     149     -0.098    -0.507
    174     150      0.488    -0.169
    175     151     -0.293     0.507
    176     152      0.195     0.338
    177     417    417       417    
    178     153     -0.293    -0.169
    179     154      0.195    -0.338
    180     155      0.488     0.169
    181     156     -0.293    -0.507
    182     157     -0.488     0.169
    183     158      0.293    -0.507
    184     159      0.098     0.507
    185     160     -0.195    -0.338
    186     161      0.585     0.000
    187     162      0.390     0.338
    188     418    418       418    
    189     163     -0.195     0.000
    190     164      0.000     0.000
    191     165      0.195     0.000
    192     508    508       508    
    193     509    509       509    
    194     166     -0.098     0.169
    195     167     -0.098    -0.169
    196     168      0.098     0.169
    197     169      0.098    -0.169
    198     419    419       419    
    199     170     -0.390     0.338
    200     171     -0.390    -0.338
    201     172      0.293     0.507
    202     173      0.293    -0.169
    203     174     -0.390     0.000
    204     175      0.098    -0.507
    205     176      0.000     0.338
    206     177      0.293     0.169
    207     178     -0.293     0.169
    208     179      0.000    -0.338
    209     420    420       420    
    210     180     -0.098     0.507
    211     181      0.390    -0.338
    212     182     -0.488    -0.169
    213     183      0.390     0.000
    214     184     -0.195     0.338
    215     185     -1.950     0.000
    216     186      0.585     1.689
    217     187     -1.170     1.351
    218     188      1.267     0.844
    219     189     -1.365     0.676
    220     421    421       421    
    221     190     -0.195     1.689
    222     191      1.170     1.351
    223     192     -1.658     0.169
    224     193      1.755     0.338
    225     510    510       510    
    226     194     -0.780     1.689
    227     195      0.293     1.520
    228     196     -1.170     1.013
    229     197      0.878     1.520
    230     422    422       422    
    231     198     -1.658     0.507
    232     199      1.462     0.844
    233     200     -1.072     1.520
    234     201     -0.293     1.520
    235     202      1.852     0.169
    236     203     -1.462     0.844
    237     204      1.462     0.507
    238     205     -1.852     0.169
    239     206      0.195     1.689
    240     207      1.072     1.182
    241     423    423       423    
    242     208      0.975     1.689
    243     209     -0.585     1.689
    244     210     -1.072     1.182
    245     211      1.658     0.169
    246     212      1.365     1.013
    247     213     -1.560     0.676
    248     214      1.560     0.676
    249     215      0.975     1.351
    250     216     -1.755     0.000
    251     217      0.000     1.689
    252     424    424       424    
    253     218     -0.780     1.351
    254     219     -1.365     0.338
    255     220      0.585     1.351
    256     221     -1.267     1.182
    257     511    511       511    
    258     222      1.170     0.676
    259     223     -0.195     1.351
    260     224      1.267     1.182
    261     225     -1.560     0.000
    262     425    425       425    
    263     226      0.780     1.689
    264     227      1.462     0.169
    265     228     -1.072     0.844
    266     229     -0.683     1.520
    267     230     -1.755     0.338
    268     231      0.975     1.013
    269     232      0.390     1.689
    270     233      1.658     0.507
    271     234      0.195     1.351
    272     235     -1.365     1.013
    273     426    426       426    
    274     236     -0.975     1.689
    275     237      1.072     1.520
    276     238     -0.390     1.689
    277     239      1.365     0.338
    278     240     -1.462     0.507
    279     241      0.683     1.520
    280     242     -0.975     1.351
    281     243      1.365     0.676
    282     244     -0.390     1.351
    283     245     -1.462     0.169
    284     427    427       427    
    285     246      0.098     1.520
    286     247     -1.170     0.676
    287     248      0.878     1.182
    288     249      1.560     0.338
    289     512    512       512    
    290     250     -0.878     1.182
    291     251      0.390     1.351
    292     252      1.072     0.844
    293     253     -1.560     0.338
    294     428    428       428    
    295     254     -0.488     1.520
    296     255      0.780     1.351
    297     256     -1.267     0.844
    298     257      0.000     1.351
    299     258      1.267     0.507
    300     259     -0.878     1.520
    301     260     -1.267     0.507
    302     261      1.170     1.013
    303     262      0.488     1.520
    304     263     -0.585     1.351
    305     429    429       429    
    306     264     -0.975     1.013
    307     265     -0.098     1.520
    308     266     -1.365     0.000
    309     267      1.072     0.507
    310     268      0.488     1.182
    311     269     -0.293     1.182
    312     270     -0.975     0.676
    313     271      0.878     0.844
    314     272     -1.072     0.169
    315     273      1.267     0.169
    316     430    430       430    
    317     274      0.195     1.013
    318     275     -0.780     1.013
    319     276      0.683     1.182
    320     277     -0.390     1.013
    321     513    513       513    
    322     278      0.878     0.507
    323     279     -1.072     0.507
    324     280      0.098     1.182
    325     281      0.488     0.844
    326     431    431       431    
    327     282     -1.170     0.000
    328     283      1.072     0.169
    329     284     -0.683     1.182
    330     285      0.975     0.676
    331     286     -0.585     0.676
    332     287     -0.098     0.844
    333     288     -0.878     0.169
    334     289      0.780     1.013
    335     290      0.780     0.338
    336     291     -0.878     0.844
    337     432    432       432    
    338     292      0.293     1.182
    339     293     -1.267     0.169
    340     294     -0.488     1.182
    341     295      0.780     0.676
    342     296     -0.878     0.507
    343     297      0.195     0.676
    344     298     -0.098     1.182
    345     299      1.170     0.338
    346     300      0.585     1.013
    347     301     -0.975     0.000
    348     433    433       433    
    349     302     -0.488     0.844
    350     303     -1.170     0.338
    351     304      0.683     0.507
    352     305     -0.585     0.338
    353     514    514       514    
    354     306      0.000     1.013
    355     307      0.683     0.844
    356     308     -0.683     0.844
    357     309      0.390     1.013
    358     434    434       434    
    359     310     -0.195     0.676
    360     311     -0.975     0.338
    361     312      0.585     0.338
    362     313     -0.780     0.676
    363     314      0.390     0.676
    364     315     -0.585     1.013
    365     316     -0.195     1.013
    366     317     -0.780     0.338
    367     318      0.098     0.844
    368     319      0.975     0.338
    369     435    435       435    
    370     320     -0.390     0.676
    371     321      0.488     0.507
    372     322     -0.488     0.507
    373     323      0.293     0.844
    374     324      0.878     0.169
    375     325     -0.780     0.000
    376     326     -0.293     0.844
    377     327      0.000     0.676
    378     328      0.585     0.676
    379     329     -0.683     0.507
    380     436    436       436    
    381     330     -0.683     0.169
    382     331      0.683     0.169
    383     999    999       999    
    384     515    515       515    
