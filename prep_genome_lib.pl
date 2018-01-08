#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/lib");
use Pipeliner;
use Cwd;
use File::Path;

my $CPU = 4;
my $outTmpDir=undef;

my $max_readlength = 150;

my $output_dir = "ctat_genome_lib_build_dir";

my $usage = <<__EOUSAGE__;

##################################################################################
#
#  Required by: STAR-Fusion and FusionInspector
#
#  --genome_fa  <string>           genome fasta file
#
#  --gtf <string>                  transcript structure annotation
#                                     Note: can restrict to coding genes and lncRNAs
#
#  --blast_pairs <string>          transcript blastn results 
#                                  in BLAST+ '-outfmt 6'  format and gzipped!
#                                     Note: gene symbols must replace transcript IDs.
#
# Required by STAR
#
#  --max_readlength <int>          max length for an individual RNA-Seq read (ie. default: $max_readlength)
#
#  Misc options:
#
#  --output_dir <string>           output directory (default: $output_dir)
#
#  --fusion_annot_lib <string>     fusion annotation library (key/val pairs, tab-delimited)
#
#  --CPU <int>                     number of threads (defalt: $CPU)
#
#  --gmap_build                    include gmap_build (for use w/ DISCASM/GMAP-fusion)
# 
#  --outTmpDir	<string>	   passed to STAR (very useful if local disks are faster than network disks)
#
##################################################################################


__EOUSAGE__

    ;



=advanced_usage

#  Experimental - when using FusionInspector with GSNAP and/or HISAT
#
#  --count_kmers                   flag to include additional build steps required by count_kmers
#                                     (requires 'jellyfish' software exist in your PATH setting)
#
#  --cdna_fa <string>              cdna fasta file
#                                     Note: header format must be:
#                                     transcript_id(tab)gene_id(tab)gene_symbol
#                                  (requires bowtie(v1)-build in PATH setting)
#


=cut

my $help_flag;
my $genome_fa_file;
my $cdna_fa_file;
my $gtf_file;
my $blast_pairs_file;
my $count_kmers;
my $fusion_annot_lib;

my $gmap_build_flag = 0;

&GetOptions ( 'h' => \$help_flag,

              # required for STAR-Fusion, FusionInspector, GMAP-fusion
              'genome_fa=s' => \$genome_fa_file,
              'gtf=s' => \$gtf_file,        
              
              # required for STAR-Fusion, FusionInspector, and GMAP-Fusion
              'blast_pairs=s' => \$blast_pairs_file,

              # required for star
              'max_readlength=i' => \$max_readlength,
              
              # optional
              'output_dir=s' => \$output_dir,
              'CPU=i' => \$CPU,
              'outTmpDir=s' => \$outTmpDir,
    
              # required for FusionInspector w/ gsnap and/or hisat
              'count_kmers' => \$count_kmers,
              'cdna_fa=s' => \$cdna_fa_file,

              # for discasm
              'gmap_build' => \$gmap_build_flag,
    
              'fusion_annot_lib=s' => \$fusion_annot_lib,
    );


if ($help_flag) {
    die $usage;
}

unless ($genome_fa_file && $gtf_file && $max_readlength && $blast_pairs_file) {
    die $usage;
}

if ($blast_pairs_file) {
    unless ($blast_pairs_file =~ /\.gz$/) {
        die "Error, file $blast_pairs_file must be gzipped-compressed";
    }
}

$genome_fa_file = &Pipeliner::ensure_full_path($genome_fa_file) if $genome_fa_file;
$cdna_fa_file = &Pipeliner::ensure_full_path($cdna_fa_file) if $cdna_fa_file;
$gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;
$blast_pairs_file = &Pipeliner::ensure_full_path($blast_pairs_file) if $blast_pairs_file;
$output_dir = &Pipeliner::ensure_full_path($output_dir) if $output_dir;
$fusion_annot_lib = &Pipeliner::ensure_full_path($fusion_annot_lib) if $fusion_annot_lib;

my $UTILDIR = $FindBin::Bin . "/util";

unless ($output_dir) {
    $output_dir = cwd();
}

my @tools_required = qw(STAR);
if ($gmap_build_flag) {
    push (@tools_required, 'gmap_build');
}
if ($cdna_fa_file) {
    push (@tools_required, 'bowtie-build');
}

my $missing_tool_flag = 0;
foreach my $tool (@tools_required) {
    my $path = `which $tool`;
    unless ($path =~ /\w/) {
        print STDERR "Error, cannot locate required tool: $tool\n";
        $missing_tool_flag = 1;
    }
}
if ($missing_tool_flag) {
    die "missing at least one required tool";
}

main: {

    my $pipeliner = new Pipeliner(-verbose => 2);
    
    #################
    # Prep the genome

    unless (-d $output_dir) {
        mkpath($output_dir) or die "Error, cannot mkpath $output_dir";
    }
    my $checkpoints_dir = "$output_dir/__checkpoints";
    unless (-d $checkpoints_dir) {
        mkpath($checkpoints_dir);
    }
    
    
    my $cmd;
    
    unless (-e "$output_dir/ref_genome.fa") {
        $cmd = "cp $genome_fa_file $output_dir/ref_genome.fa";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_ref_genome.fa.ok"));
    }

    
    ###############################
    ## and copy the annotation file
    
    unless (-e "$output_dir/ref_annot.gtf") {

        $cmd = "cp $gtf_file $output_dir/ref_annot.gtf";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_ref_annot.gtf.ok"));
    }
    
    unless (-e "$output_dir/ref_annot.gtf.mini.sortu") {
        $cmd = "bash -c \" set -eof pipefail; $UTILDIR/gtf_to_exon_gene_records.pl $output_dir/ref_annot.gtf  | sort -k 1,1 -k4,4g -k5,5g | uniq  > $output_dir/ref_annot.gtf.mini.sortu \" ";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_sort_ref_annot_gtf_mini.ok"));
        
    }
    
        
    # build star index
    my $star_index = "$output_dir/ref_genome.fa.star.idx";
    unless (-d $star_index) {
        mkpath $star_index or die "Error, cannot mkdir $star_index";
    }

    my $maybe_tmpdir= defined($outTmpDir)? " --outTmpDir $outTmpDir " : "";

    $cmd = "STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index $maybe_tmpdir "
            . " --genomeFastaFiles $output_dir/ref_genome.fa "
            . " --limitGenomeGenerateRAM 40419136213 "
            . " --genomeChrBinNbits 16 " # needed for >4k contigs w/ FI
            . " --sjdbGTFfile $gtf_file "
            . " --sjdbOverhang $max_readlength ";
    
    $pipeliner->add_commands(new Command($cmd, "$star_index/build.ok"));

    unless (-e "$output_dir/ref_annot.gtf.gene_spans") {
        $cmd = "$UTILDIR/gtf_to_gene_spans.pl $output_dir/ref_annot.gtf > $output_dir/ref_annot.gtf.gene_spans";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_ref_annot.gtf.gene_spans.ok") );
    }
    

    #######################
    # index the blast pairs
    
    if ($blast_pairs_file) {
        $cmd = "$UTILDIR/index_blast_pairs.pl $output_dir/blast_pairs.idx $blast_pairs_file ";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_blast_pairs.idx.ok"));

        # remove blast pairs between genes that physically overlap on the genome
        $cmd = "$UTILDIR/index_blast_pairs.remove_overlapping_genes.pl $output_dir";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_prune_blast_pairs_overlapping_genes.ok"));
    }
    

    ####################################
    ## integrate protein structure info
    
    $cmd = "$UTILDIR/build_prot_info_db.pl --gtf $gtf_file --genome_fa $genome_fa_file --out_prefix $output_dir/ref_annot";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_prot_info_db.ok"));
    
    
    ######################
    # build the fusion annotation database
    # 
    $cmd = "$UTILDIR/build_fusion_annot_db_index.pl --gene_spans $output_dir/ref_annot.gtf.gene_spans --out_db_file $output_dir/fusion_annot_lib.idx";
    if ($fusion_annot_lib) {
        $cmd .= " --key_pairs $fusion_annot_lib";
    }
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_fusion_annot_lib.idx.ok"));
    

    ############################################################
    ## Sections below are optional depending on optional params
    ############################################################
    
    if ($gmap_build_flag) {
        
        #########################
        # build GMAP genome index
        
        $cmd = "gmap_build -D $output_dir -d ref_genome.fa.gmap -k 13 $output_dir/ref_genome.fa";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_ref_genome.fa.gmap.ok"));
    }
    
    if ($cdna_fa_file) {

        
        ## This is for FusionInspector, in case HISAT or GSNAP is used.
        

        ##########################
        # Prep the cdna fasta file
        
        unless (-e "$output_dir/ref_cdna.fasta") {
            $cmd = "cp $cdna_fa_file $output_dir/ref_cdna.fasta";
            $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_ref_cdna.fasta.ok"));
        }
                
        # index the fasta file
        $cmd = "$UTILDIR/index_cdna_seqs.pl $output_dir/ref_cdna.fasta";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_ref_cdna.fasta.idx.ok"));
        
        # build the bowtie index:
        $cmd = "bowtie-build $output_dir/ref_cdna.fasta $output_dir/ref_cdna.fasta";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_ref_cdna.fasta.bowtie_idx.ok"));
        
    }
                
    if ($count_kmers) {
        my $cmd = "jellyfish count -t $CPU -m 25 -s 1000000000 --canonical $output_dir/ref_genome.fa";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_jf.count.ok"));

        $cmd = "jellyfish dump -L 2 mer_counts.jf > $output_dir/ref_genome.jf.min2.kmers";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_jf.dump.ok"));

        $cmd = "rm mer_counts.jf";
        $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/_jf.cleanup.ok"));
        
    }
    
    $pipeliner->run();

    exit(0);
}




        
    
