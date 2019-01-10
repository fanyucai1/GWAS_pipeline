#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::long;
use Cwd;
use Config::IniFiles;
use FindBin qw($Bin);
use File::Basename;
my $tassel="/share/work/biosoft/tassel/tasseladmin-tassel-v5/run_pipeline.pl";
my $perl="/share/work/biosoft/perl/perl-5.22.1/bin/perl";
my($hmp,$kinship,$structure,$phenotype,$outdir);
GetOptions(
    "hmp:s"=>\$hmp,       
    "p:s"=>\$phenotype,
    "k:s"=>\$kinship,
    "c:s"=>\$structure,
    "o:s"=>\$outdir,
    );

sub usage{
    print qq {
This script use tassel to run GWAS.
usage:
#use MLM
perl $0 -hmp ##hmp.txt -p traits.txt -k kinship.txt -c structure.txt -o /path/to/diretory
options:
-hmp            the genotype data(force)
-k              the kinship data(force)
-c              the strycture data
-p              the traits data(force)
-o              output diretory
Email:fanyucai1\@126.com
2016.6.29
verion 1.0
    };
    exit;
}

#GWAS use mlm
system "echo '$perl $tassel  -fork1 -h $hmp -filterAlign -filterAlignMinFreq 0.05 -fork2 -r $phenotype -fork3 -p $structure -excludeLastTrait -fork4 -k $kinship -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm ­-mlmOutputFile $outdir/tassel -runfork1 -runfork2 -runfork3 -runfork4'>$outdir/GWAS.sh";
#tree_build
system "echo '$perl $tassel -h $hmp -tree Neighbor -treeSaveDistance false -export $outdir/tree.nj.txt '>>$outdir/GWAS.sh";

system "sh $outdir/GWAS.sh";