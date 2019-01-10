#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Cwd;
use Getopt::Long;

my ($kinship,$genotype,$phenotype,$structure,$outdir,$repeat,$FDR,$MAF);
my $R="/usr/local/bin/Rscript";
$repeat||=100;
$FDR||=1;
$MAF||=0;
GetOptions(
    "KI:s"=>\$kinship,
    "G:s"=>\$genotype,
    "P:s"=>\$phenotype,
    "CV:s"=>\$structure,
    "o:s"=>\$outdir,
    "rel:s"=>\$repeat,
    "FDR:s"=>\$FDR,
    "MAF:s"=>\$MAF,
           );

sub usage{
    print qq{
This script will run the GWAS analysis using GAPIT.
usage:
perl $0 -KI kinship.txt -G genotype.txt -P phenotype -CV structure.txt -o /path/to/directory
options:
-G                  genotypic data in either standard HapMap format(force)
-P                  phenotypes in the text file of phenotypic data(force)
-rel                Program runing time - repetition times(default:100)
-FDR                Threshold to Filter SNP on FDR(default:1)
-MAF                Minor Allele Frequency to Filter SNPs in GWAS Reports(default:0)
-o                  the output directory
    };
    exit;
}

if(!$genotype || !$outdir || !$phenotype)
{
    &usage();
}
system "mkdir -p $outdir";
system "echo '#!$R
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler)
library(\"scatterplot3d\")
source(\"http://www.zzlab.net/GAPIT/emma.txt\")
source(\"http://www.zzlab.net/GAPIT/gapit_functions.txt\")
myG<-read.table(\"$genotype\",head=FALSE)
myY<-read.table(\"$phenotype\",head=TRUE)
myGAPIT <- GAPIT(Y=myY,G=myG,
kinship.cluster=c(\"average\", \"complete\", \"ward\"),
kinship.group=c(\"Mean\", \"Max\"),
SNP.MAF=$MAF,
SNP.FDR=$FDR,
PCA.total=3,
Model.selection = TRUE)'>$outdir/gapit.Rscript";
system "cd $outdir && $R $outdir/gapit.Rscript";







