use strict;
use warnings;

open(net,"network.txt") or die "cannot open file $!";
my(@fields, @gnams, @columns,@fields2,%ghash1,%ghash2);



open(anno,"anno.txt") or die "can not open file $!";
open(OUT,">newnetwor.txt") or die "can not open $!";
while(<anno>){
chomp();
     s/\r\n?/\n/;  # remove return
     chomp();
     @fields2 = split(/\t/, $_);
	 #print "$fields2[0]\t$fields2[5]\n";
	 if($fields2[5] eq "---"){
	 $ghash2{$fields2[0]}=$fields2[0];
	 }else{
	 $ghash2{$fields2[0]}=$fields2[5];}

}



while(<net>){
 chomp();
     s/\r\n?/\n/;  # remove return
     chomp();
     @fields = split(/\t/, $_);
	 #print "$fields[1]\t$fields[2]\n";
	 if(exists $ghash2{$fields[1]} or exists $ghash2{$fields[2]}){
	 print OUT "$fields[0]\t$fields[1]\t$ghash2{$fields[1]}\t$fields[2]\t$ghash2{$fields[2]}\t$fields[3]\n";
	 
	 }else{print "$_\n";}
	 
}

close(net);


