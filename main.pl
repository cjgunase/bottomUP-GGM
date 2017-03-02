use lib '/home/cjgunase/myModules/share/perl5/';#Stat TRY
use lib '/home/cjgunase/myModules/lib64/perl5/';#PDLLIST

use Statistics::Multtest qw(:all);
use Try::Tiny;
use PDL;
use Statistics::RankCorrelation;
use Statistics::Distributions;
use List::MoreUtils qw/ uniq /;
use Benchmark qw(:hireswallclock);

#---------time spent------------
my $starttime = Benchmark->new;
my $finishtime;
my $timespent;
open(TIME,">time.txt") or die;
#---------------------------------

majic("At_stem_rma_july2012_AGI22PW_final.csv","selected_rowsum_top_300_TFs.csv","layer1",22);
qx(sort -nrk 4,4 temp1.txt>temp.txt);
network("selected_rowsum_top_300_TFs.csv","file0genes1.csv","file0genes2.csv","layer1.txt",25);

majic("file0genes1.csv","file0genes2.csv","layer2",22);
qx(sort -nrk 4,4 temp1.txt>temp.txt);
network("file0genes2.csv","file1genes1.csv","files1genes2.csv","layer2.txt",20);

majic("file1genes1.csv","files1genes2.csv","layer3",22);
qx(sort -nrk 4,4 temp1.txt>temp.txt);
network("files1genes2.csv","file2genes1.csv","file2genes2.csv","layer3.txt",15);
#
majic("file2genes1.csv","file2genes2.csv","layer4",22);
qx(sort -nrk 4,4 temp1.txt>temp.txt);
network("file2genes2.csv","file3genes1.csv","file3genes2.csv","layer4.txt",10);

majic("file3genes1.csv","file3genes2.csv","layer5",22);
qx(sort -nrk 4,4 temp1.txt>temp.txt);
network("file3genes2.csv","file4genes1.csv","file4genes2.csv","layer5.txt",5);





sub network{
	my($file,$genefile1,$genefile2,$layer,$number)=@_;
	
	open(LAY,"temp.txt") or die "temp File not found $!";
 	open(OUT,">$layer") or die "not found";
  	my $layer1=$number;
  	my($count,%layer1hash);
  	#<LAY>;
  	
  	while(<LAY>){

my @values = split('\t', $_);
$layer1hash{$values[2]}="";
$count = keys %layer1hash;
if($count eq $layer1){ last;}

}
close(LAY);

open(LAY,"temp.txt") or die "temp File not found $!"; 
#<LAY>;
while(<LAY>){
	s/\r\n?/\n/;  # remove return
	chomp();
	@values = split('\t', $_);
	if (exists $layer1hash{$values[2]}){ 
		print OUT "$_\n";
	}
	
}
close(LAY);
open(TF,"$file") or die "File 1 not found $!";
open(file1,">$genefile1") or die "file not found";
open(file2,">$genefile2") or die "file not found";


while(<TF>){
	s/\r\n?/\n/;  # remove return
	chomp();
my ($k, $v)=split(",",$_,2);

if(exists $layer1hash{$k}){
	
	print file1 "$_\n";
}else{
	print file2 "$_\n";
}

}




}
sub majic{

	my($file1,$file2,$whatlayer,$n)=@_;
	
my($count,$count1,%hoa);
open(F1,$file1) or die "File 1 not found $!";#f1=pw gene file
open(F2,$file2) or die "File 2 not found $!";#f2=TF file
open(OUT,">temp1.txt") or die "file not found $!";
#open(OUT2,">temp2.txt") or die "file not found $!";
#print OUT2 "PWgene\tTF1\tTF2\tcorrx_1\tp_valuex_1\tdiff1\tp_value1\tcorrx_2\tp_valuex_2\tdiff2\tp_value\tinterferece\n";
						
print OUT "layer\tPWgene\tTF\tFreq\n";
	
	
	while(<F1>){
	$count++;
	my ($k, $v)=split(",",$_,2);
	$f1gene{$count}=$k;
	$f1value{$count}=$v;
	}
	
	while(<F2>){
	$count1++;
	my ($k, $v)=split(",",$_,2);
	$f2gene{$count1}=$k;
	$f2value{$count1}=$v;
	}
	
	foreach (sort keys %f1gene) {
	$hoa{$f1gene{$_}};
	my @arr_varx=$f1value{$_};
	my @varx=split(",",@arr_varx[0]);
	my $vx=\@varx;
		
		for($a=1;$a<$count1;$a++){					#a<TFlength
		
			for($b=$a+1;$b<($count1+1);$b++){				#b<TFlenght+1
			
			my @arr_var1=$f2value{$a};
			my @arr_var2=$f2value{$b};
			@var1=split(",",@arr_var1[0]);
			@var2=split(",",@arr_var2[0]);
			$v1=\@var1;
			$v2=\@var2;
			
			my $corrx_1= corr_func($vx,$v1);
			my $corrx_2= corr_func($vx,$v2);
			
			my $tstatx_1= tstat($corrx_1,$n,2);
			my $p_valuex_1=Statistics::Distributions::tprob ($n-2,$tstatx_1);#one tailed put df and t statistics to calculate the P-values FOR COR(X,Y)
		 	$p_valuex_1=$p_valuex_1*2;# p_value = for correlation X1
		 	
		 	my $tstatx_2= tstat($corrx_2,$n,2);
			my $p_valuex_2=Statistics::Distributions::tprob ($n-2,$tstatx_2);#one tailed put df and t statistics to calculate the P-values FOR COR(X,Y)
		 	$p_valuex_2=$p_valuex_2*2;# p_value = for correlation X2
		 	
		 	my ($diff1,$diff2,$diffse1,$diffse2,$ISEXCEP)=spc_func($vx,$v1,$v2,$n);
		 	if($ISEXCEP > 0){ next;}
			my $z1=$diff1/$diffse1;
			my $p_value1=Statistics::Distributions::uprob ($z1);#one tailed put df and t statistics to calculate the P-values
		 	#$p_value1=$p_value1*2;#COR(1,3)-CORX1,(3,2)
		 	my $z2=$diff2/$diffse2;
			my $p_value2=Statistics::Distributions::uprob ($z2);#one tailed put df and t statistics to calculate the P-values
		 	#$p_value2=$p_value2*2;#COR(2,3)-CORX2,(3,1)
				
				if($p_valuex_1<0.05){
					if($p_value1<0.05){
						$c1++; #correlation and difference both significant
						push @{ $hoa{$f1gene{$_}} }, "$f2gene{$a}", "$f2gene{$b}";
						#print OUT2 "$f1gene{$_}\t$f2gene{$a}\t$f2gene{$b}\t$corrx_1\t$p_valuex_1\t$diff1\t$p_value1\t$corrx_2\t$p_valuex_2\t$diff2\t$p_value2\tBoth_interfere";
						
						}else{
										$c2++;#correlation significant difference not significance
										
									}
																		
								}else{
									
									if ($p_valuex_2 < 0.05){
										
										if($p_value2 < 0.05){$c3++; #other correlation significant difference significant
											push @{ $hoa{$f1gene{$_}} }, "$f2gene{$a}", "$f2gene{$b}";
											#print OUT2 "$f1gene{$_}\t$f2gene{$a}\t$f2gene{$b}\t$corrx_1\t$p_valuex_1\t$diff1\t$p_value1\t$corrx_2\t$p_valuex_2\t$diff2\t$p_value2\tBoth_interfere";
										}else{$c4++; #correlation significant difference not siginificant
										}
										
									}else{
										$c5++; #correlation not siginicant difference not signicant
										#print OUT2 "$f1gene{$_}\t$f2gene{$a}\t$f2gene{$b}\t$corrx_1\t$p_valuex_1\t$diff1\t$p_value1\t$corrx_2\t$p_valuex_2\t$diff2\t$p_value2\tBoth_NOT_interfere";
									}
									
								}
								
								}
	
	
							}	
					}
					#print OUT2 "X12cor_diff_sig\tcorr_sig_diff_not\tX21cor_diff_sig	\ttcorr_sig_diff_not	\tboth_not_significant	\n";
					#print OUT2 "$c1\t			$c2\t	$c3\t	$c4\t	$c5\n";
					
		foreach $gene ( keys %hoa ) {
			my @array = @{ $hoa{$gene} };
    #print "$gene:". scalar @array."\n";
    my %freq;
$freq{$_}++ for @array;
    
    foreach (sort keys %freq) {
    print OUT "$whatlayer\t$gene\t$_\t$freq{$_}\n";
  }
    
    }				
	
	
}
sub corr_func {
	    my ($v1,$v2) = @_;
        try {
        my $cor1_2 = Statistics::RankCorrelation->new( $v1, $v2, sorted => 1 );
  			 $cor1_2 = $cor1_2->spearman;
        return ($cor1_2);
        
			}catch{        
       				return (0);
        		}
}
sub tstat{
	my ($v,$n,$df)=@_;
	try{
	$tstat=($v*sqrt($n-$df))/sqrt(1-$v*$v);
	return ($tstat);
}catch{
	return 0;
}
}
sub spc_func {
	    my ($v1,$v2,$v3,$n) = @_;
        try {
        	
         my $cor1_2 = Statistics::RankCorrelation->new( $v1, $v2, sorted => 1 );
  			 $cor1_2 = $cor1_2->spearman;
  			
  			my $cor3_1 = Statistics::RankCorrelation->new( $v1, $v3, sorted => 1 );
  			 $cor3_1 = $cor3_1->spearman;
  			
  			my $cor3_2 = Statistics::RankCorrelation->new( $v3, $v2, sorted => 1 );
  			 $cor3_2 = $cor3_2->spearman;
  			
         
         my $cor1_32=($cor3_1-($cor1_2*$cor3_2))/sqrt(1-($cor3_2*$cor3_2));
         
         my $cor2_31=($cor3_2-($cor3_1*$cor1_2))/sqrt(1-($cor3_1*$cor3_1));
         
         my $diff1=abs($cor3_1-$cor1_32);
         my $diff2=abs($cor3_2-$cor2_31);
         
         
         my $var_r_12=((1-$cor1_2**2)**2)/$n;
         my $var_r_23=((1-$cor3_2**2)**2)/$n;
         my $var_r_13=((1-$cor3_1**2)**2)/$n;
         
         my $covar13_12=(0.5*(2*$cor3_2-$cor3_1*$cor1_2)*(1-$cor3_2**2-$cor3_1**2-$cor1_2**2)+$cor3_2**3)/$n;
         my $covar13_32=(0.5*(2*$cor1_2-$cor3_1*$cor3_2)*(1-$cor3_2**2-$cor3_1**2-$cor1_2**2)+$cor1_2**3)/$n;
         my $covar23_12=(0.5*(2*$cor3_1-$cor3_2*$cor1_2)*(1-$cor3_2**2-$cor3_1**2-$cor1_2**2)+$cor3_1**3)/$n;
         
         
         my $pd1=$cor3_2/sqrt(1-($cor3_2*$cor3_2));#r1_2
         my $pd2=($cor1_2+($cor3_1*$cor3_2))/sqrt(1-($cor3_2*$cor3_2));#r2_3
         my $pd3=1-(1/sqrt($cor3_2*$cor3_2));#r3_1
         
         my $a = pdl [
             		[ $pd1, $pd2, $pd3 ]
    					];
    	
    	 my $b = pdl [
    				[ $var_r_12, $covar23_12, $covar13_12 ],
    				[ $covar23_12, $var_r_23, $covar13_32 ],
    				[ $covar13_12, $covar13_32, $var_r_13 ]
						];
						
		 my $c = pdl [
             		[ $pd1 ], 
             		[ $pd2 ], 
             		[ $pd3 ] 
    					];
    									
         
         my $var_diff1=$a x $b x $c;
         $var_diff1=at($var_diff1,0,0);
	     my $diffse1=sqrt($var_diff1);
         
         
         $pd1=$cor3_1/sqrt(1-($cor3_1*$cor3_1));#r1_2
         $pd2=($cor1_2+($cor3_1*$cor3_2))/sqrt(1-($cor3_1*$cor3_1));#r2_3
         $pd3=1-(1/sqrt($cor3_1*$cor3_1));#r3_1
         
         $a = pdl [
             		[ $pd1, $pd2, $pd3 ]
    					];
    	
    	 $b = pdl [
    				[ $var_r_12, $covar23_12, $covar13_12 ],
    				[ $covar23_12, $var_r_23, $covar13_32 ],
    				[ $covar13_12, $covar13_32, $var_r_13 ]
						];
						
		 $c = pdl [
             		[ $pd1 ], 
             		[ $pd2 ], 
             		[ $pd3 ] 
    					];
    									
         
         my $var_diff2=$a x $b x $c;
         $var_diff2=at($var_diff2,0,0);
	     my $diffse2=sqrt($var_diff2);
         
         
         return ($diff1,$diff2,$diffse1,$diffse2);
        
         #return ($cor1_2,$cor3_1,$cor3_2,$cor1_32,$cor2_31)
} catch {
        #print "exception\n";
       return (0,0,1,1,1,1);
        
}
         
         
    }
#--------------------------------------------------

#---------------------Calculate the time spent---------------------
	$finishtime = Benchmark->new;
	$timespent = timediff($finishtime,$starttime);
	print TIME "\nDone!\nSpent ". timestr($timespent);