#! /usr/bin/perl
use warnings;
use strict;
use Cwd;
use Config::IniFiles;
use FindBin qw($Bin);
use File::Basename;
use GetOpt::long;

my $R="/share/work/biosoft/R/R-v3.0.0/bin/";
my($kinship,$genotype_hmp,$phenotype,$structure,$outdir);
GetOptions(
    "G:s"=>\$genotype_hmp,
    "Y:s"=>\$phenotype,
    "o:s"=>\$outdir,
           );

sub usage{
    print qq{
This script use gapit to GWAS.
usage:
perl $0 -K kinship.txt -G hmp.txt -Y phenotype.txt -Q structure.txt -o /path/to/directory
options:
-G           Genotype Data in Hapmap Format
-Y           phenotype
-o           output directory
Email:fanyucai1\@126.com
2016.6.28
    };
    exit;
}

system "echo'
#!$R/Rscript
ibrary(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler)
source(\"http://zzlab.net/GAPIT/gapit_functions.txt\")
source(\"http://www.zzlab.net/GAPIT/emma.txt\")'>$outdir/GWAS.sh";

if ($kinship)
{
    system "echo 'myKI <- read.table(\"$kinship\", head = FALSE)'>>$outdir/GWAS.sh";
}

if ($phenotype)
{
    system "echo 'myY <- read.table(\"$phenotype\", head = TRUE)'>>$outdir/GWAS.sh";
}

if ($genotype_hmp)
{
    system "echo 'myG <- read.table(\"$genotype_hmp\" , head = FALSE)'>>$outdir/GWAS.sh";
}

if ($structure)
{
    system "echo 'myCV <- read.table(\"$structure\", head = TRUE)'>>$outdir/GWAS.sh";
}






