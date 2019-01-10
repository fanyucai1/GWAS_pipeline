#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;

my $R="/share/work/biosoft/R/R-v3.0.0/bin/";
my $plink="/share/work/biosoft/plink/plink_v1.9/plink";
my($input,$outdir);
my $window||=500;
my $MAF||=0.05;
my $geno||=0.1;
my $hwe||=0.0001;
my $snp ||=1000000;
my $r ||=0.05;
my $step||=100;
GetOptions(
    "i:s"=>\$input,
    "o:s"=>\$outdir,
    "w:s"=>\$window,
    "MAF:s"=>\$MAF,
    "geno:s"=>\$geno,
    "hwe:s"=>\$hwe,
    "r:s"=>\$r,
    "bin:s"=>\$step,
           );

sub usage{
    print qq{
This script caculate the LD and Plot Linkage Disequilibrium Decay.
usage:
perl $0 -i plink -o /path/to/directory
options:
-i                  the file of prefix input plink (http://pngu.mgh.harvard.edu/~purcell/plink/ld.shtml)
-w                  Average r2 was calculated for pairwise markers in a window and averaged across the whole genome. (default 1000kb),or set 500(kb)
-o                  the output directory(force)
-MAF                Minimum minor allele frequency to include a marker(default:0.05)
-hwe                Exclude markers with a HW p-value smaller than <threshold>(default:0.0001),or set(0)
-geno               Exclude markers with less than <threshold> valid data.(default:0.1)
-snp                To only analyse SNPs that are not more than 10 SNPs apart, for example, use the option (default is 10 SNPs),or set 100
-r                  To report only values above a particular value (this only applies when the --r2 and not the --r command is used) (default is 0.05)
-bin                caculate the average of LD and plot point (default is 100),when genome is large you could increase this value
Email:fanyucai1\@126.com
2016.7.1
};
    exit;
}
if (!$input || !$outdir)
{
    &usage();
}
`$plink --file $input --make-founders --geno $geno --maf $MAF --r2 --ld-window-r2 $r --ld-window $snp --ld-window-kb $window --out $outdir/test_plink`;
`awk \'{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}\' $outdir/test_plink.ld >$outdir/ld_matrix`;
open(IN,"$outdir/ld_matrix");
open(OUT,">$outdir/bin.txt");
my %hash1;
my %hash2;

while (<IN>)
{
    chomp;
    if ($_!~/CHR/)
    {
        my @array=split(/ /,$_);
        $hash1{int(($array[4]-$array[1])/$step)}+=$array[6];
        $hash2{int(($array[4]-$array[1])/$step)}++; 
    } 
}
print OUT "position\tR2\n";
foreach  my $key(sort{$a<=>$b} keys %hash1)
{
    $hash1{$key}=$hash1{$key}/$hash2{$key};
    print OUT $key,"\t",$hash1{$key},"\n";
}

system "echo '
#!$R/Rscript
library(ggplot2)
a<-read.table(\"$outdir/bin.txt\",sep=\"\\t\",header=T)
png(\"$outdir/ld_decay.png\",res=300,width=3000,height=3000)
ggplot(a,aes(position,R2))+geom_point(colour=\"#990000\")+scale_y_continuous(expression(r^2),breaks = c(0,0.2,0.4,0.6,0.8,1),limits = c(0, 1))+theme(axis.text= element_text(size=16),axis.title= element_text(size=20))+scale_x_continuous(\"Pairwise distance(Kb)\",breaks = c(0,100000,200000,300000,400000,500000), labels = c(\"0\",\"100\", \"200\", \"300\",\"400\",\"500\"))
dev.off()
pdf(\"$outdir/ld_decay.pdf\",width=10,height=8)
ggplot(a,aes(position,R2))+geom_point(colour=\"#990000\")+scale_y_continuous(expression(r^2),breaks = c(0,0.2,0.4,0.6,0.8,1),limits = c(0, 1))+theme(axis.text= element_text(size=16),axis.title= element_text(size=20))+scale_x_continuous(\"Pairwise distance(Kb)\",breaks = c(0,100000,200000,300000,400000,500000), labels = c(\"0\",\"200\", \"400\", \"600\",\"800\",\"1000\"))
dev.off()
'>$outdir/plink_ld.Rscript";
`$R/Rscript $outdir/plink_ld.Rscript`;
