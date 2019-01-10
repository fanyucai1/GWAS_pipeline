#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my $Annovar="/share/work/biosoft/annovar2016_06_24/annovar/";
my $gtfToGenePred="/share/work/fanyc/pipeline/GWAS_pipeline/sub/gtfToGenePred";
my $perl="/share/work/biosoft/perl/perl-5.22.1/bin/perl";
my $qsub="/usr/local/bin/qsub-sge.pl";
my($gtf,$gff,$outdir,$vcf,$ref,$prefix);

GetOptions(
    "gff:s"=>\$gff,
    "gtf:s"=>\$gtf,
    "vcf:s"=>\$vcf,
    "o:s"=>\$outdir,
    "ref:s"=>\$ref,
    "prefix:s"=>\$prefix,
           );

sub usage{
    print qq{
This script was used to annotation vcf file by Annovar.
usage:
perl $0 -gtf ###.gtf -vcf ###.vcf -o /path/to/outdir -ref ref.fa -prefix ##
options:
-gtf                input gff (gtf) file(force)
-vcf                input vcf file (force)
-o                  output directory(force)
-ref                the reference sequence(force)
-prefix             the prefix of output(force)
Email:fanyucai1\@126.com
2016.6.24
    };
    exit;
}
if (!$vcf || !$outdir || !$ref ||!$prefix)
{
    &usage();
}

`mkdir -p $outdir`;
system "echo 'cd $outdir/'>$outdir/annovar.sh";
if ($gtf)
{
    system "echo '$gtfToGenePred -genePredExt $gtf $prefix\_refGene.txt'>>$outdir/annovar.sh";
    `ln -s $gtf $outdir/`;
}
`ln -s $ref $outdir/`;
system "echo '$perl $Annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile $ref $prefix\_refGene.txt --out $prefix\_refGeneMrna.fa'>>$outdir/annovar.sh";
system "echo '$perl $Annovar/table_annovar.pl $vcf $outdir/ --vcfinput --out final --buildver $prefix --remove --protocol refGene --operation g'>>$outdir/annovar.sh";
my $line=`wc -l $outdir/annovar.sh`;
chomp($line);
if (`hostname` !~ /cluster/)
{
    system "sh $outdir/annovar.sh";
}
else
{
    `perl $qsub -line $line`;
}
