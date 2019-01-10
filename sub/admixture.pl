#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;
my($bed,$outdir,$kcluster);
my $R="/share/work/biosoft/R/R-v3.2.3/bin/";
my $admixture="/share/work/biosoft/admixture/admixture_linux-1.3.0/admixture";
$outdir||=getcwd;
GetOptions(
    "bed:s"=>\$bed,
    "k:s"=>\$kcluster,
    "o:s"=>\$outdir,
        );

sub usage{
    print qq{
This script will find the best cluster and plot stucture by admixture.
usage:
perl $0 -bed plink.bed -k 10 -o /path/to/outdir
options:
-bed                the output from plink will be same directory with plink.bim and plink.fam files
-k                  the process will find best number from 1 to kcluster
-xlab               the name of xlab
-ylab               the name of ylab
-o                  the output directory
    };
    exit;
}

if (!$bed || !$outdir)
{
    &usage();
}


my @cluster;
my $basename=basename($bed,  ".bed");
system "echo 'k\terro'>$outdir/k_erro.txt";
my $pwd=`pwd`;
chomp($pwd);
for(my $k=1;$k<=$kcluster;$k++)
{
    system "$admixture --cv $bed $k >$outdir/log_$k.out";
    $cluster[$k-1]=`grep CV $outdir/log_$k.out |awk -F\: \'{print \$2}\'`;
    chomp($cluster[$k-1]);
    system "echo '$k\t$cluster[$k-1]'>>$outdir/k_erro.txt";
}
my @sorted = sort { $a <=> $b } @cluster;
print $sorted[0];
for(my $k=0;$k<$kcluster;$k++)
{
    if ($cluster[$k]==$sorted[0])
    {
        my $best=$k+1;
        system "echo '
        #!$R/Rscript
        tbl=read.table(\"$outdir/$basename.$best.Q\")
        png(filename=\"$outdir/admixture.png\",height=1000,width=3000)
        barplot(t(as.matrix(tbl)), col=rainbow($best), xlab=\"Individual\",ylab=\"Ancestry\", border=NA)
        dev.off()
        pdf(file=\"$outdir/admixture.pdf\")
        barplot(t(as.matrix(tbl)), col=rainbow($best), xlab=\"Individual\",ylab=\"Ancestry\", border=NA)
        dev.off()
        a=read.table(\"$outdir/k_erro.txt\",header=T)
        library(ggplot2)
        png(filename=\"$outdir/erro.png\",height=1000,width=1500)
        ggplot(a,aes(k,erro))+geom_point()+ geom_smooth()
        dev.off()
        pdf(file=\"$outdir/erro.pdf\")
        ggplot(a,aes(k,erro))+geom_point()+ geom_smooth()
        dev.off()
        '>$outdir/admixture.R";
    }   
}
system "$R/Rscript $outdir/admixture.R";
