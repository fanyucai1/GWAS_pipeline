#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $R="/share/work/biosoft/R/R-v3.2.3/bin/";
my ($input,$outdir);
GetOptions(
    "i:s"=>\$input,
    "o:s"=>\$outdir,  
           );

sub usage
{
    print qq{
This script plot manhattan and qq plot.
usage:
perl $0 -i input.txt -o /path/to/directory
-i          the input file,this file looks like this:
                SNP CHR BP         P
                rs1   1  1 0.9148060
                rs2   1  2 0.9370754
                rs3   1  3 0.2861395
                rs4   1  4 0.8304476
                rs5   1  5 0.6417455
                rs6   1  6 0.5190959
-o          the output directory
    };
    exit;
}
if (!$input || !$outdir)
{
    &usage();
}


system "echo '
#!$R/Rscript
library(qqman)
a<-read.table(\"$input\",header=T,sep=\"\\t\")
as.dataframe(a)
png(\"$outdir/manhattan.png\",res=300,height=3000,width=3000)
manhattan(a, col = c(\"blue4\", \"orange3\"))
dev.off()
pdf(\"$outdir/manhattan.pdf\")
manhattan(a, col = c(\"blue4\", \"orange3\"))
dev.off()
png(\"$outdir\QQ.png\",res=300,height=3000,width=3000)
qq(a\$P, main = \"Q-Q plot of GWAS p-values\")
dev.off()
pdf(\"$outdir\QQ.pdf\")
qq(a\$P, main = \"Q-Q plot of GWAS p-values\")
dev.off()
'>>$outdir/manhattan_qq.Rscript";

system "$R/Rscript $outdir/manhattan_qq.Rscript";

