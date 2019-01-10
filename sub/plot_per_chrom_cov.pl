#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $samtools="/share/work/biosoft/samtools/samtools-1.3/samtools";
my $R="/share/work/biosoft/R/R-v3.2.3/bin/Rscript";
my $bedtools="/share/work/biosoft/bedtools/bedtools2-v2.25.0/bin/bedtools";
my ($bam,$genomeSize,$outdir,$prefix);
my $window||=10000;
GetOptions(
    "bam:s"=>\$bam,
    "genomeSize:s"=>\$genomeSize,
    "w:s"=>\$window,
    "p:s"=>\$prefix,
    "o:s"=>\$outdir,
           );

sub usage{
    print qq{
This script plot the coverage per chromosome.
options:
-bam                input the sort bam(force)
-genomeSize         input bed file contains chromosome size this file from samtools faidx ouput(force)
-w                  the silde window(default:10000)
-p                  the prefix of output(force)
-o                  the output directory(force)
usage:
perl $0 -bam sample.sort.bam -genomeSize genome.fa.fai -o /path/to/outdir/ -p sample
Email:fanyc\@biomics.com.cn
2016.3.15
    };
    exit;
}
if (!$bam || !$genomeSize|| !$outdir || !$prefix) {
    &usage();
}

system "$bedtools makewindows -g $genomeSize -w $window >$outdir/$prefix.$window.bed";
open(IN,"$outdir/$prefix.$window.bed");
open(OUT,">$outdir/$prefix.$window\_chr.bed");
while (<IN>)
{
    chomp;
    if ($_!~"Scaffold")
    {
        print OUT "$_\n";
    }
}

system "$samtools bedcov $outdir/$prefix.$window\_chr.bed $bam >$outdir/$prefix.$window.cov";
system "echo '#!$R'>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'x<-read.table(\"$outdir/$prefix.$window.cov\")'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'chr<-x[,1]'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'pos<-x[,3]'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'depth<-x[,4]/$window'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'depth<-log2(depth)'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'df  <- data.frame(chr, pos, depth)'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'require(ggplot2)'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'pdf(file=\"$outdir/$prefix\_cov.pdf\",width=40,height=30)'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'p <- ggplot(data = df, aes(x=pos, y=depth),binwidth = 0.1) + geom_area(aes(fill=chr))'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'p + facet_wrap(~ chr, ncol=1)'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'dev.off()'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'png(filename=\"$outdir/$prefix\_cov.png\",width=6000,res=300,height=5000)'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'p <- ggplot(data = df, aes(x=pos, y=depth),binwidth = 0.1) + geom_area(aes(fill=chr))'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'p + facet_wrap(~ chr, ncol=1)'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "echo 'dev.off()'>>$outdir/$prefix.plot_per_chrom_cov.R";
system "$R $outdir/$prefix.plot_per_chrom_cov.R"
