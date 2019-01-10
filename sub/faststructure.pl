#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename;

my $python="/share/work/biosoft/python/Python-v2.7.11/bin/python";
my $faststructure="/share/work/biosoft/faststructer/fastStructure/";
my $R="/share/work/biosoft/R/R-v3.2.3/bin/Rscript";
my $qsub="/usr/local/bin/qsub-sge.pl";
my $gcta="/share/work/biosoft/gcta/gcta64";
my($bed,$outdir,$prefix,$pdfw,$pdfh,$pngh,$pngw,$k,$dir);
$pdfh||=5;
$pdfw||=10;
$pngw||=3000;
$pngh||=1500;
$k||=10;
GetOptions(
    "bed:s"=>\$bed,
    "o:s"=>\$outdir,
    "p:s"=>\$prefix,
    "pdfw:s"=>\$pdfw,
    "pdfh:s"=>\$pdfh,
    "pngw:s"=>\$pngw,
    "pngh:s"=>\$pngh,
    "dir:s"=>\$dir,
    "k:s"=>\$k,
           );

sub usage{
    print qq{
This script will given a sturcture of community.
usage:
perl $0 -bed plink -o /path/to/dir/ -p prefix -dir /path/to/diectory/
options:
-bed            the prefix of bed file,plink.bed, plink.bim and plink.fam files are all present in the same path(force)
-dir            the directory of bed file(force)
-o              the output directory(force)
-p              the prefix of output(force)
-pdfh           the height of pdf(default:8)
-pdfq           the width of pdf(default:12)
-pngh           the height of png(default:1500)
-pngw           the width of png(default:3000)
-k              the max cluster in community(default:10)
Email:fanyucai1\@126.com
2016.8.30
    };
    exit;
}
if(!$bed || !$outdir || !$prefix || !$k || !$dir)
{
    &usage();
}
system "mkdir -p $outdir";
my $dirname=dirname($bed);
my $basename=basename($bed);
open(CLU,">$outdir/faststructure.sh");
for(my $j=2;$j<=$k;$j++)
{
    print CLU "cd $dir && $python $faststructure/structure.py --input=$bed -K $j --output=$outdir/$prefix --full --seed=100\n";
}
if(`hostname`=~"cluster")
{
    system "perl $qsub $outdir/faststructure.sh";
}
else
{
    &usage();
}
`awk '{print \$2}' $bed.fam >$outdir/$prefix.sampleID.txt`;
my (@sample,$num);
$num=0;
open(SAM,"$outdir/$prefix.sampleID.txt");
while(<SAM>)
{
    chomp;
    $sample[$num]=$_;
    $num++; 
}

for(my $j=1;$j<=$k;$j++)
{
    open(IN,"$outdir/$prefix.$j.meanQ");
    open(OUT,">$outdir/$prefix.$j.cluster");
    print OUT "sample\tcluster\tvalue\n";
    my $num=0;
    while(<IN>)
    {
        chomp;
        my @array=split(/  /,$_);
        for(my $k=1;$k<=$j;$k++)
        {
            my $number=$array[$k-1]*100;
            print OUT "$sample[$num]\tcluster$k\t$number\n";
        }
        $num++;
    }
    close IN;
    close OUT;
}

`echo "#!$R">$outdir/cluster.Rscript`;

for(my $j=1;$j<=$k;$j++)
{
    system"echo '
x<-read.table(\"$outdir/$prefix.$j.cluster\",header=T,sep=\"\\t\")
library(ggplot2)
p=ggplot(x,aes(sample,value,fill=cluster))+ geom_bar(stat=\"identity\")+theme(axis.text.x  = element_text(angle=30, vjust=0.5, size=7))+ylab(\"Percent(%)\")+xlab(\"sampleID\")
png(\"$outdir/$prefix.$j.png\",res=300,height=$pngh,width=$pngw)
p
dev.off()
pdf(\"$outdir/$prefix.$j.pdf\",height=$pdfh,width=$pdfw)
p
dev.off()
    '>>$outdir/cluster.Rscript";
}

system "$python $faststructure/chooseK.py --input=$outdir/$prefix >$outdir/best_cluster.txt";
=head
open(IN,"$outdir/best_cluster.txt");
my $numt;
while(<IN>)
{
    chomp;
    my @array=split(/=/,$_);
    $numt=$array[1];
}
system "$gcta --bfile $bed --make-grm --autosome --out $bed && $gcta --grm $bed --pca 3 --out $outdir/$num.pcatmp";
=cut
system "$R $outdir/cluster.Rscript";