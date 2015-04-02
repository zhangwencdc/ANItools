#!/usr/bin/perl

=head1 Name

=head1 Description

=head1 

  Author: zhangwen, zhangwen@icdc.cn
  Version: 1.0, Date: 2015-02-01

=head1 Usage:

perl  ANI_one_to_matrix.pl -I genome.fa -O ./ -Genus Streptococcus -Tag 05ZYH33
or
perl  ANI_one_to_matrix.pl -I genome.fa -O ./ -Species Streptococcus_suis -Tag 05ZYH33

=head1 Example

=cut
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use warnings;
use Pod::Text;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;
##get options from command line into variables and set default values
###-----------------Time record-----------------------------###
my $Time_Start = sub_format_datetime(localtime(time())); #......
my $Data_Vision = substr($Time_Start,0,10);
my ($genome,$outdir,$genus,$species,$tag,$HELP);
GetOptions(
		"I:s"=>\$genome,
		"O:s"=>\$outdir,
		"Genus:s"=>\$genus,
		"Species:s"=>\$species,
		"Tag:s"=>\$tag,
		"help"=>\$HELP
);
die `pod2text $0` if ($HELP || !defined $genome || !defined $tag);
die `pod2text $0` if (!defined $genus && !defined $species);
if(defined $species){$species=~tr/_/ /;}
my %config;
parse_config("$Bin/config.txt",\%config);
my $glimmer=$config{glimmer};
my $glimmer_csh=$config{glimmer_csh};
my $program=$config{ANItools};
my $formatdb=$config{formatdb};#Blast建库路径
my $blast=$config{blast};
my $trex=$config{trex};
my $qsub=$program."/qsub-sge.pl";
my $blast_analysis=$program."/blast_analysis.pl";
my $iden=0.6;
my $len2=0.7;
##Glimmer 基因预测
my $time=sub_format_datetime(localtime(time()));
print "$time\n";
print "Gene prediction\n";
system "mkdir $outdir/glimmer\n";

#print "sh $Bin/glimmer_train.sh $glimmer linear $genome $outdir/glimmer/$tag $glimmer_csh\n";
system "sh $Bin/glimmer_train.sh $glimmer linear $genome $outdir/glimmer/$tag $glimmer_csh\n";
#print "$glimmer  -o50 -g110 -t30 -l $genome $outdir/glimmer/$tag.icm  $outdir/glimmer/$tag\n";
system "$glimmer/glimmer3  -o50 -g110 -t30 -l $genome $outdir/glimmer/$tag.icm  $outdir/glimmer/$tag\n";
#print "perl $Bin/predict_convert.pl --predict glimmer  --finalname $tag --log --verbose $outdir/glimmer/$tag.predict $genome\n";
system "perl $Bin/predict_convert.pl --predict glimmer  --finalname $tag --log --verbose $outdir/glimmer/$tag.predict $genome\n";

$time=sub_format_datetime(localtime(time()));
print "$time\n";
print "Blast\n";
###读取genome list明确target
my $list=$Bin."/all_bac.list"  ;#genome list
open (LIST,"$list")||die "Can't open FILE:$list\n";
my %target;
while(1){
	my $line=<LIST>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	if((defined $genus && $a[0] eq $genus)||(defined $species && $a[1] eq $species)){
		my $stest=basename($a[3]);
		my $sname=$a[2].",".substr($stest,0,length($stest)-4);
		$target{$sname}=$a[3];
	}
}
my @target=keys %target;
close LIST;
###计算ANI
my $cds=$outdir."/glimmer/".$tag.".predict.cds";
my $strain1=$tag;my $s1=$strain1;
my $query=$cds;
my %ANI;
my $blast_data=$outdir;
my %gene;
foreach my $strain2 (@target) {
		my $target=$target{$strain2};
		my $s2=$strain2;
		if($query eq "" || $query eq "NA"||$target eq "" || $target eq "NA"){$ANI{$strain1}{$strain2}="NA";next;}
		###blast
		$b=$outdir."/".$s1.".".$s2.".blast";
		my $parser=$b.".parser";
		my $bdata=$target.".nhr";
		my $target_name=basename($target);
		unless(-e $bdata){
		system "$formatdb -i $target -p F -o T -n $target_name\n ";
		system "$blast -p blastn -d $target_name -i $query -o $outdir/$s1.$s2.blast -X 150 -q -1 -F F\; ";
		}else{
			system "$formatdb -i $target -p F -o T -n $target\n ";
		system "$blast -p blastn -d $target -i $query -o $outdir/$s1.$s2.blast -X 150 -q -1 -F F\; ";
		}
		system "perl $program/blast_parser.pl $b >$parser\n";
		
		#######
		####blast parser

		open(FILE,"$parser")||die "Can't open File: $parser\n";
		
		my $query_len;
		my $align_len;
		my $identity;
		my %query_len;
		my %align_len;
		my %identity;
		while (<FILE>) {
				chomp;
				my $line=$_;
				my @a=split"\t",$line;
				if($a[0]=~/Query/){next;}
				my $q_len=$a[1];
				my $i=$a[8];
				my $p=$a[9];
				my $align=$a[11];
				my $q_anno=$a[14];
				my $t_anno=$a[15];
				$q_len=~s/\,//;
				$align=~s/\,//;
				#unless($identity=~/\d/ && $positive=~/\d/){print "$identity,$positive,$align,$q_len,$line\n"};
				unless($q_len>0){next;}
						#print "$i,$align,$q_len,$line\n";
				if(defined $iden){unless($i>=$iden){next;}}
				#print "$i,$align,$q_len,$line\n";
				if(defined $len2){unless(($align/$q_len)>=$len2){next;}}

				if(exists $identity{$a[0]} and $identity{$a[0]}>($align*$i)){next;} #排除重复序列
				$gene{$a[0]}++;
				$query_len{$a[0]}=$q_len;
				$align_len{$a[0]}=$align;
				$identity{$a[0]}=$align*$i;
				
		}
		close FILE;
#		system "rm $parser\n";
		my @key=sort keys %query_len;
		foreach my $key (@key) {
				$query_len+=$query_len{$key};
				$align_len+=$align_len{$key};
				$identity+=$identity{$key};
				
		}
		 ###排除染色体和质粒比对的情况
		 #print "$align_len,$query_len,$out1\n";
		unless($query_len>0){next;}
		 unless(($align_len/$query_len)>=0.5){ next;}
		 my $align=$align_len/$query_len;
		 my $i=$identity/$align_len;
		$ANI{$strain1}{$strain2}=$i;
		#print "$strain1,$strain2,$align_len,$query_len,$identity,$i\n";

}
##读取原有的matrix列表内容
my $matrix=$program."/ANI_matrix";
my $target_genus;
if(defined $genus){$target_genus=$genus;}elsif(defined $species){my @species=split" ",$species;$target_genus=$species[0];}
my $ANI_matrix=$matrix."/".$target_genus.".ANI";
open(FILE,"$ANI_matrix")||die "Can't open FILE:$ANI_matrix\n";
my @name;
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	my $num=@a;
	if($a[0] eq "ANI"){
		foreach (1..($num-1)) {
			$name[$_]=$a[$_];
		}
	}else{
		my $name2=$a[0];
		foreach  (1..($num-1)) {
			$ANI{$name2}{$name[$_]}=$a[$_];
		}
	}

}
close FILE;
$ANI{$tag}{$tag}=1;
###结果输出ANI matrix
push @target,$tag;
open(OUT,">$outdir/$tag.ANI");
print OUT "ANI\t";
foreach my $strain (@target) {
		my @strain1=split",",$strain;
		my $strain_n1;
		if($strain1[0] eq $tag){$strain_n1=$strain1[0];}else{$strain_n1=$strain1[0]."__".$strain1[1];}
		print OUT "$strain_n1\t";

}
print OUT "\n";
foreach my $strain1 (@target) {

	my @strain1=split",",$strain1;
	my $strain_n1;
	if($strain1[0] eq $tag){$strain_n1=$strain1[0];}else{$strain_n1=$strain1[0]."__".$strain1[1];}
		
		my $aa_value="";
		my $aa;
		foreach my $strain2 (@target) {
			my $ANI_value;
		#	my @strain2=split",",$strain2;
		#	$strain2=$strain2[0];
			if($strain2 eq $tag){$ANI_value=$ANI{$strain2}{$strain1};
			}else{
			$ANI_value=$ANI{$strain1}{$strain2};
			}
			$aa_value=$aa_value.$ANI_value."\t";
			if($ANI_value=~/[0-9a-zA-Z]/){$aa++;}
		}
		if($aa>1){  ###含有至少一个有效数据
		print OUT "$strain_n1\t";
		print OUT "$aa_value";
		print OUT "\n";
		}
}
close OUT,
###ANI matrix 转换为Distance matrix
$time=sub_format_datetime(localtime(time()));
print "$time\n";
print "Tree1\n";
system "perl $program/ANItoDistance.pl $outdir/$tag.ANI $outdir/$tag.distance\n";
##distance matrix转换为tree 生成*.nwk和mytree.svg文件
system "perl $program/distancetotree.pl $outdir/$tag.distance $outdir/$tag.nwk\n";

###输出query中特异的gene 若identity>=60% align>70%，则认为有比对基因
$time=sub_format_datetime(localtime(time()));
print "$time\n";
print "Run Speical Gene\n";
open(FILE,"$cds")||die "Can't open FILE:$cds\n";
open(OUT,">$outdir/$tag.specialgene");
my $mark=0;
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	if(substr($line,0,1) eq ">"){
	my @a=split" ",$line;
	my $name=substr($a[0],1);
	if(exists $gene{$name}){$mark=0;}else{$mark=1;}
	}
	if($mark==1){print OUT "$line\n";}
}
close FILE;
close OUT;

####trex #第二种生成tree的方法
$time=sub_format_datetime(localtime(time()));
print "$time\n";
print "Run Tree2\n";
system "perl $program/ANItoDistance2.pl $outdir/$tag.ANI $outdir/$tag.2.distance\n";
system "$trex/trex -method=2 -option=1 -k=5 -p=0 -option1=1 -optionFunction=1 -newickFile=$outdir/$tag.tree2.nwk  -outmatrix=$outdir/matrixfile.txt -dataType=1 -bootstrap=0 -stat=output.txt -optionTR=3 -input=$outdir/$tag.2.distance -output=treeFile.txt -modele=$trex/data/modele.txt \n";
print "Finished\n";
my $Time_End= sub_format_datetime(localtime(time()));
print "Running from [$Time_Start] to [$Time_End]\n";



#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+



sub sub_format_datetime #.....
{
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon, $day, $hour, $min, $sec);
}
#####
sub gc{
	my $seq=shift @_;
	my $gc=$seq=~tr/(G|C|g|c)/(G|C|g|c)/;
	my $l=length($seq);


	return ($gc,$l);
}

##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		next if(/#/);
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;

		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}
