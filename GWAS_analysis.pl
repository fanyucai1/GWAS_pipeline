#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Config::IniFiles;
use FindBin qw($Bin);
use File::Basename;
my ($input);
my $thread ||=5;
my $minPruning ||=2;
my $MAF||=0.05;
my $geno||=0.1;
my $cluster||=10;
my $prune ||=0.3;
my $window||=10000;
my $coverage||=5;
my $mbp||=20;
my $hwe||=0.000001;
my $max_DP||=1000;
my $min_DP ||=2;
my $seq||="all.q";
my $ass_maxproc||="10";
my $del;
my $GATK_thread||=3;
GetOptions(
    "i:s"=>\$input,
    "t:s"=>\$thread,
    "w:s"=>\$window,
    "minPruning:s"=>\$minPruning,
    "c:s"=>\$coverage,
    "MAF:s"=>\$MAF,
    "geno:s"=>\$geno,
    "r:s"=>\$prune,
    "k:s"=>\$cluster,
    "mbp:s"=>\$mbp,
    "del:s"=>\$del,
        );

if (!$input)
{
    &usage();
}
sub usage
{
    print qq(
This script will run the gwas pipeline.
example:
perl $0 -i config.ini
Options:
-i                  the config file
-o                  the output directory
-del                excludes sites with indentifiers matching,only give coverage per chromosome <forexample: scaffold>
-w                  the silde window plot coverage per chromosome(default:10000)
-t                  number of threads,default:5
-minPruning         Minimum support to not prune paths in the graph(GATK HaplotypeCaller Advanced Parameters,default:2)
-geno               Missing rates(default:0.05)
-MAF                Minor Allele Frequency(default:0.05,or 0.03,0.01)
-mbp                Minimum base quality required to consider a base for calling(default:20)
-hwe                Removing SNPs out of Hardy-Weinberg equilibrium(default:0.000001)
-r                  Linkage disequilibrium based SNP pruning(default:0.3 or 0.5,http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune)
-k                  the cluster of k used by admixture(default:10)
    );
    exit;   
}
####################################qsub
sub qsub()
{
	my ($shfile, $queue, $ass_maxproc) = @_ ;
    $queue||="all.q";
    $ass_maxproc||=10;
    my $cmd = "perl /usr/local/bin/qsub-sge.pl --maxproc $ass_maxproc --queue $queue --resource vf=15G --reqsub $shfile --independent" ;
	if (`hostname` !~ /cluster/)
    {
		print "Only cluster could run this process ";
        exit;
	}
    my $flag=system($cmd);
    if($flag !=0)
    {
		die "qsub [$shfile] die with error : $cmd \n";
        exit;
	}
}
sub run_or_die()
{
	my ($cmd) = @_ ;
	my $flag = system($cmd) ;
	if ($flag != 0){
		print "die with: $cmd\n";
		exit;
	}
}
####################################read the configure and make output_directory
my %ini;
tie %ini, 'Config::IniFiles', ( -file => $input);
foreach my $key (keys %{$ini{dir}})
{
    $ini{dir}{$key}=~ s/\$outdir/$ini{dir}{outdir}/;
    system "mkdir -p $ini{dir}{$key}";
}
######################################mapping
my $strings="";
my $count=0;
my $bamS;
foreach my $key (keys %{$ini{"rawdata"}})
{
    my @array=split(/\|/,$ini{"rawdata"}{$key});
    #align raw data use bwa
    if ($count==0)
    {
        system "echo '$ini{software}{bwa} mem -M -t $thread  $ini{database}{ref} $array[0] $array[1] > $ini{dir}{Mapping}/$key-pe.sam'>$ini{dir}{shell}/bwa_align_1.sh";
        #sam2bam
        system "echo '$ini{software}{samtools} view -bS  $ini{dir}{Mapping}/$key-pe.sam -o $ini{dir}{Mapping}/$key-pe.bam'>$ini{dir}{shell}/bwa_align_2.sh";
        #Add read groups, coordinate sort and index using AddOrReplaceReadGroups
        system "echo 'rm $ini{dir}{Mapping}/$key-pe.sam && $ini{software}{java} -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{picard} AddOrReplaceReadGroups INPUT=$ini{dir}{Mapping}/$key-pe.bam OUTPUT=$ini{dir}{Mapping}/$key\_addRG.bam RGID=$key RGLB=$key RGPL=illumina RGPU=$key RGSM=$key SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR=/share/work/tmp'>$ini{dir}{shell}/bwa_align_3.sh";
        #marking duplication
        system "echo 'rm $ini{dir}{Mapping}/$key-pe.bam && $ini{software}{java} -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{picard} MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=$ini{dir}{Mapping}/$key\_addRG.bam OUTPUT=$ini{dir}{Mapping}/$key\_markdup.bam METRICS_FILE=$ini{dir}{Mapping}/$key\_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=/share/work/tmp'>$ini{dir}{shell}/bwa_align_4.sh";
        #Create a target list of intervals to be realigned
        system "echo '$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T RealignerTargetCreator -allowPotentiallyMisencodedQuals -fixMisencodedQuals -R $ini{database}{ref} -I $ini{dir}{Mapping}/$key\_markdup.bam -o $ini{dir}{Mapping}/$key\_realignment_targets.list'>$ini{dir}{shell}/bwa_align_5.sh";
        #Perform realignment of the target intervals
        system "echo '$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T IndelRealigner -R $ini{database}{ref} -I $ini{dir}{Mapping}/$key\_markdup.bam -allowPotentiallyMisencodedQuals -targetIntervals $ini{dir}{Mapping}/$key\_realignment_targets.list -o $ini{dir}{Mapping}/$key\_realigned_reads.bam'>$ini{dir}{shell}/bwa_align_6.sh";
        #Variant-only calling on DNAseq(https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)
        system "echo '$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -R $ini{database}{ref} -T HaplotypeCaller -nct $GATK_thread -I $ini{dir}{Mapping}/$key\_realigned_reads.bam -mbq $mbp --minPruning $minPruning --emitRefConfidence GVCF -o $ini{dir}{Variant}/$key.g.vcf'>$ini{dir}{shell}/bwa_align_7.sh";
    }
    else
    {
        system "echo '$ini{software}{bwa} mem -M -t $thread $ini{database}{ref} $array[0] $array[1] > $ini{dir}{Mapping}/$key-pe.sam'>>$ini{dir}{shell}/bwa_align_1.sh";
        #sam2bam
        system "echo '$ini{software}{samtools} view -bS  $ini{dir}{Mapping}/$key-pe.sam -o $ini{dir}{Mapping}/$key-pe.bam'>>$ini{dir}{shell}/bwa_align_2.sh";
        #Add read groups, coordinate sort and index using AddOrReplaceReadGroups
        system "echo 'rm $ini{dir}{Mapping}/$key-pe.sam && $ini{software}{java} -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{picard} AddOrReplaceReadGroups INPUT=$ini{dir}{Mapping}/$key-pe.bam OUTPUT=$ini{dir}{Mapping}/$key\_addRG.bam RGID=$key RGLB=$key RGPL=illumina RGPU=$key RGSM=$key SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR=/share/work/tmp'>>$ini{dir}{shell}/bwa_align_3.sh";
        #marking duplication
        system "echo 'rm $ini{dir}{Mapping}/$key-pe.bam && $ini{software}{java} -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{picard} MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=$ini{dir}{Mapping}/$key\_addRG.bam OUTPUT=$ini{dir}{Mapping}/$key\_markdup.bam METRICS_FILE=$ini{dir}{Mapping}/$key\_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=/share/work/tmp'>>$ini{dir}{shell}/bwa_align_4.sh";
        #Create a target list of intervals to be realigned
        system "echo '$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T RealignerTargetCreator -allowPotentiallyMisencodedQuals -fixMisencodedQuals -R $ini{database}{ref} -I $ini{dir}{Mapping}/$key\_markdup.bam -o $ini{dir}{Mapping}/$key\_realignment_targets.list'>>$ini{dir}{shell}/bwa_align_5.sh";
        #Perform realignment of the target intervals
        system "echo '$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T IndelRealigner -R $ini{database}{ref} -I $ini{dir}{Mapping}/$key\_markdup.bam -allowPotentiallyMisencodedQuals -targetIntervals $ini{dir}{Mapping}/$key\_realignment_targets.list -o $ini{dir}{Mapping}/$key\_realigned_reads.bam'>>$ini{dir}{shell}/bwa_align_6.sh";
        #Variant-only calling on DNAseq(https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)
        system "echo '$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx10g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -R $ini{database}{ref} -T HaplotypeCaller -nct $GATK_thread -I $ini{dir}{Mapping}/$key\_realigned_reads.bam -mbq $mbp --minPruning $minPruning --emitRefConfidence GVCF -o $ini{dir}{Variant}/$key.g.vcf'>>$ini{dir}{shell}/bwa_align_7.sh";
    }
    $strings.="--variant $ini{dir}{Variant}/$key.g.vcf ";
    $count++;
    $bamS.="$ini{dir}{Mapping}/$key\_markdup.bam ";
}
system "echo 'perl $ini{software}{BamS} -bam $bamS -o $ini{dir}{Mapping}/ -del $del -genomeSize $ini{database}{ref_fai} -w $window'>$ini{dir}{shell}/bwa_statistics.sh";
##############################
    &qsub("$ini{dir}{shell}/bwa_align_1.sh");
    &qsub("$ini{dir}{shell}/bwa_align_2.sh");
    &qsub("$ini{dir}{shell}/bwa_align_3.sh");
    &qsub("$ini{dir}{shell}/bwa_align_4.sh");
    &qsub("$ini{dir}{shell}/bwa_align_5.sh");
    &qsub("$ini{dir}{shell}/bwa_align_6.sh");
    &qsub("$ini{dir}{shell}/bwa_align_7.sh");
    &qsub("$ini{dir}{shell}/bwa_statistics8.sh");
###############################

#Perform joint genotyping on gVCF files produced by HaplotypeCaller
    
    open(VCF,"$ini{dir}{shell}/get_VCF_9.sh");
    print VCF "$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T GenotypeGVCFs -nct $GATK_thread -R $ini{database}{ref} $strings -o $ini{dir}{Variant}/snp_indel.vcf";
    &qsub("$ini{dir}{shell}/get_VCF.sh","great.q");
    
#Extract the SNPs from the call set and #Apply the filter to the SNP call set(http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)
    open(SNP,">$ini{dir}{shell}/$ini{dir}{shell}/get_SNP_10.sh");
    print SNP "$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T SelectVariants -nct $GATK_thread -R $ini{database}{ref} -V $ini{dir}{Variant}/snp_indel.vcf -selectType SNP -o $ini{dir}{Variant}/raw_snps.vcf && ";
    print SNP "$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T VariantFiltration -nct $GATK_thread -R $ini{database}{ref} -o $ini{dir}{Variant}/filtered_snps.vcf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --clusterWindowSize 5 --clusterSize 2 --filterName \"my_snp_filter\" --variant $ini{dir}{Variant}/raw_snps.vcf && ";
    print SNP "$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T SelectVariants -nct $GATK_thread -R $ini{database}{ref} -V $ini{dir}{Variant}/filtered_snps.vcf -o $ini{dir}{Variant}/snps.vcf -ef";
    &qsub("$ini{dir}{shell}/get_SNP_10.sh","great.q");
    
#Extract the Indels from the call set and Apply the filter to the Indel call set
    open(INDEL,">$ini{dir}{shell}/get_indel_11.sh");
    print INDEL "$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T SelectVariants -nct $GATK_thread -R $ini{database}{ref} -V $ini{dir}{Variant}/snp_indel.vcf -selectType INDEL --maxIndelSize 30 -o $ini{dir}{Variant}/raw_indels.vcf && ";
    print INDEL "$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T VariantFiltration -nct $GATK_thread -R $ini{database}{ref} -o $ini{dir}{Variant}/filtered_indels.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"my_indel_filter\" -V $ini{dir}{Variant}/raw_indels.vcf && ";
    print INDEL "$ini{software}{java} -Djava.io.tmpdir=/share/work/tmp -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{gatk} -T SelectVariants -nct $GATK_thread -R $ini{database}{ref} -V $ini{dir}{Variant}/filtered_indels.vcf -o $ini{dir}{Variant}/indels.vcf -ef";
    &qsub("$ini{dir}{shell}/get_VCF_11.sh","great.q");
###################################

#imputation of ungenotyped markers
    open(IMPUTE,">$ini{dir}{shell}/impute_12.sh");
    print IMPUTE "cd $ini{dir}{imputation} && $ini{software}{java} -Xmx30g  -XX:ParallelGCThreads=$thread -jar $ini{software}{Beagle} gtgl=$ini{dir}{Variant}/snps.vcf out=impute.gtgl";
    &qsub("$$ini{dir}{shell}/impute_12.sh","great.q");

#filter_vcf
    open(VCF2,">$ini{dir}{shell}/filter_vcf_13.sh");
    print VCF2 "cd $ini{dir}{file_format} && $ini{software}{vcftools} --gzvcf $ini{dir}{imputation}/impute.gtgl.vcf.gz --maf $MAF --max-missing $geno --recode --out plink";
    &qsub("$ini{dir}{shell}/filter_vcf_13.sh");

#build nj tree
    open(TREE,">$ini{dir}{shell}/tree_14.sh");
    print TREE "cd $ini{dir}{file_format} && $ini{software}{vk} phylo tree nj $ini{dir}{file_format}/plink.recode.vcf 1>nj.tre 2>nj.log";
    &qsub("$ini{dir}{shell}/tree_14.sh");
    
#vcf2ped
    open(PED,">$ini{dir}{shell}/vcf2plink_15.sh");
    print PED "cd $ini{dir}{plink} && $ini{software}{vcftools} --vcf $ini{dir}{file_format}/plink.recode.vcf --plink --out plink";
    &qsub("$ini{dir}{shell}/vcf2plink_15.sh");
    
#ped2bed
    open(BED,">$ini{dir}{shell}/ped2bed_16.sh");
    print BED "cd $ini{dir}{plink} && $ini{software}{plink} --noweb --make-bed --file plink --out plink";
    &qsub("$ini{dir}{shell}/ped2bed_16.sh");
    
#bed2structure using admixture
    open(STRU,">$ini{dir}{shell}/admixture_17.sh");
    print STRU "cd $ini{dir}{plink} && $ini{software}{perl} $Bin/sub/admixture.pl -bed $ini{dir}{Variant}/plink.bed -k $cluster -o $ini{dir}{admixture}";
    &qsub("$ini{dir}{shell}/admixture_17.sh");
    
#ped2Hapmap(https://groups.google.com/forum/#!topic/tassel/Ra8vYilaKmI)
    open(HAP,">$ini{dir}{shell}/plink2tassel_18.sh");
    print HAP "cd $ini{dir}{plink} && $ini{software}{perl} $ini{software}{tassel}/run_pipeline.pl -importGuess plink.ped -export plink -exportType Hapmap";
    &qsub("$ini{dir}{shell}/plink2tassel_18.sh");

#Linkage disequilibrium based SNP pruning(http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune)
    open(PRUNE,">$ini{dir}{shell}/prune_19.sh");
    print PRUNE "cd $ini{dir}{plink} &&  $ini{software}{plink} --noweb --bfile $ini{dir}{Variant}/plink --hwe $hwe --indep-pairwise 100 10 0.3 && ";
    print PRUNE "cd $ini{dir}{plink} &&  $ini{software}{plink} --noweb --bfile $ini{dir}{Variant}/plink --extract $ini{dir}{Variant}/plink.prune.in --make-bed --out $ini{dir}{Variant}/prune";
    &qsub("$ini{dir}{shell}/prune_19.sh");
    
#use annovar to annotate the snp
    open(ANN,">$ini{dir}{shell}/annovar_20.sh");
    print ANN "$ini{software}{perl} $ini{software}{annovar} -gtf $ini{database}{gtf} -vcf $ini{dir}{Variant}/snps.vcf -o $ini{dir}{annovar}/SNP -ref $ini{database}{ref} -prefix SNP_Anno\n";
    print ANN "$ini{software}{perl} $ini{software}{annovar} -gtf $ini{database}{gtf} -vcf $ini{dir}{Variant}/indels.vcf -o $ini{dir}{annovar}/indel -ref $ini{database}{ref} -prefix indel_Anno";
    &qsub("$ini{dir}{shell}/annovar_20.sh");

#use tassel to given kinship matrix(https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Kinship/Kinship)
    open(KIN,">$ini{dir}{shell}/kinship_21.sh");
    print KIN "$ini{software}{perl} $ini{software}{tassel}/run_pipeline.pl -importGuess $ini{dir}{file_format}/QC.hmp.txt -KinshipPlugin -method Dominance_Normalized_IBS -endPlugin -export $ini{dir}{kinship}/kinship.txt -exportType SqrMatrix\n";
    &qsub("$ini{dir}{shell}/kinship_21.sh");
 
#plink2LD(input ped file)
    open(LD,">$ini{dir}{shell}/LD_22.sh");
    print LD "cd $ini{dir}{plink} && perl $Bin/sub/plink_LD.pl -i plink -o $ini{dir}{LD}";
    &qsub("$ini{dir}{shell}/LD_22.sh");
    
#faststructure
    open(FAST,">$ini{dir}{shell}/faststructure_23.sh");
    print FAST "perl $Bin/sub/faststructure.pl -bed plink -o $ini{dir}{fast} -p fast -dir $ini{dir}{plink}";
    &qsub("$ini{dir}{shell}/faststructure_23.sh");
    
#PCA(pcatmp.eigenvec)
    open(PCA,">$ini{dir}{shell}/PCA_24.sh");
    print PCA "cd $ini{dir}{plink} && $ini{software}{gcta64} --bfile plink --make-grm --autosome --out temp && ";
    print PCA "$ini{software}{gcta64} --grm temp --pca 3 --out pca";
    &qsub("$ini{dir}{shell}/PCA_24.sh");
    