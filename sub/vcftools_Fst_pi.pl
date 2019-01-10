#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my $vcftools="/share/work/biosoft/vcftools/vcftools_v0.1.14/bin/";
my $env="export PERL5LIB=/share/work/biosoft/vcftools/vcftools_v0.1.14/share/perl5/:\$PERL5LIB";
my ($vcf,$outdir,$pop1,$pop2);
my $win ||=100000;
my $step ||=10000;
GetOptions(
    "i:s"=>\$vcf,
    "o:s"=>\$outdir,
    "win:s"=>\$win,       
    "step:s"=>\$step,
    "p1:s"=>\$pop1,
    "p2:s"=>\$pop2,
           );
sub usage{
    print qq{
This script could calculate Fst and pi.
usage:
perl $0 -i input.vcf -win 100000 -step 2000 -p1 pop1.txt -p2 pop2.txt -o /path/to/outdir
-i                  the input vcf file
-win                the slide window(defualt:100K)
-step               the step(default:10K)
-p1                 The user must supply text files that contain lists of individuals (one per line) that are members of each population.
-p2                 The user must supply text files that contain lists of individuals (one per line) that are members of each population.
-o                  output directory
Email:fanyucai1\@126.com
2016.5.29
    };
    exit;
}

if (!$win || !$vcf || !$step || !$pop1 || !$outdir)
{
    &usage();
}

system "echo '$env'>$outdir/Fst.sh";
system "echo 'cd $outdir'>>$outdir/Fst.sh";
system "echo 'nohup  $vcf/vcftools  --vcf $vcf --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-step $step --fst-window-size $win --ouput $pop1.$pop2.$win.$step.Fst&'>>$outdir/Fst.sh";
system "echo 'nohup  $vcf/vcftools  --vcf $vcf --weir-fst-pop $pop1 --weir-fst-pop $pop2 --ouput $pop1.$pop2.site.Fst&'>>$outdir/Fst.sh";
system "echo 'nohup  $vcf/vcftools  --vcf $vcf --keep $pop1 --site-pi --ouput $pop1.site.pi&'>>$outdir/Fst.sh";
system "echo 'nohup  $vcf/vcftools  --vcf $vcf --keep $pop2 --site-pi --ouput $pop2.site.pi&'>>$outdir/Fst.sh";
system "echo 'nohup  $vcf/vcftools  --vcf $vcf --keep $pop1 --site-pi --window-pi $win --window-pi-step $step --ouput $pop1.$win.$step.pi&'>>$outdir/Fst.sh";
system "echo 'nohup  $vcf/vcftools  --vcf $vcf --keep $pop2 --site-pi --window-pi $win --window-pi-step $step --ouput $pop2.$win.$step.pi&'>>$outdir/Fst.sh";
