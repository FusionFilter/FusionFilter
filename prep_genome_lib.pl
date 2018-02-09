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



my $help_flag;
my $genome_fa_file;
my $gtf_file;
my $fusion_annot_lib;

my $gmap_build_flag = 0;

&GetOptions ( 'h' => \$help_flag,

              # required for STAR-Fusion, FusionInspector, GMAP-fusion
              'genome_fa=s' => \$genome_fa_file,
              'gtf=s' => \$gtf_file,        
            
              # required for star
              'max_readlength=i' => \$max_readlength,
              
              # optional
              'output_dir=s' => \$output_dir,
              'CPU=i' => \$CPU,
              'outTmpDir=s' => \$outTmpDir,
    
              # for discasm
              'gmap_build' => \$gmap_build_flag,
    
              'fusion_annot_lib=s' => \$fusion_annot_lib,
    );


if ($help_flag) {
    die $usage;
}

unless ($genome_fa_file && $gtf_file && $max_readlength) {
    die $usage;
}

my @required_tools = ("STAR", "makeblastdb", "blastn");
if ($gmap_build_flag) {
    push (@required_tools, "gmap_build");
}

&check_for_required_tools(@required_tools);


$genome_fa_file = &Pipeliner::ensure_full_path($genome_fa_file) if $genome_fa_file;
$gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;
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

    my $checkpoints_dir = "$output_dir/__chkpts";
    unless (-d $checkpoints_dir) {
        mkpath($checkpoints_dir) or die "Error, cannot mkpath $checkpoints_dir";
    }
            
    my $cmd = "cp $genome_fa_file $output_dir/ref_genome.fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir/ref_genome.fa.ok"));
        
    ###############################
    ## and copy the annotation file
    
    $cmd = "cp $gtf_file $output_dir/ref_annot.gtf";
    $pipeliner->add_commands(new Command($cmd, "$output_dir/ref_annot.gtf.ok"));
    
    # extract exon records    
    $cmd = "bash -c \" set -eof pipefail; $UTILDIR/gtf_to_exon_gene_records.pl $output_dir/ref_annot.gtf  | sort -k 1,1 -k4,4g -k5,5g | uniq  > $output_dir/ref_annot.gtf.mini.sortu \" ";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/ref_annot.gtf.mini.sortu.ok"));
            
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
    
    $cmd = "$UTILDIR/gtf_to_gene_spans.pl $output_dir/ref_annot.gtf > $output_dir/ref_annot.gtf.gene_spans";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/ref_annot.gtf.gene_spans.ok") );
    
    
    #############################
    ## Begin the BLAST compendium

    # CDS and ncRNA blastn for quick homology identification
    
    $cmd = "$UTILDIR/gtf_file_to_feature_seqs.pl $output_dir/ref_annot.gtf $output_dir/ref_genome.fa CDSplus > ref_annot.cdsplus.fa";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdsplus.fa.ok"));

    $cmd = "makeblastdb -in ref_annot.cdsplus.fa -dbtype nucl";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdsplus.fa.blidx.ok"));

    $cmd = "blastn -query ref_annot.cdsplus.fa -db ref_annot.cdsplus.fa -max_target_seqs 10000 -outfmt 6 -evalue 1e-10 -num_threads $CPU -dust no > ref_annot.cdsplus.allvsall.outfmt6";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdsplus.allvsall.outfmt6.ok"));

    
    $cmd = "bash -c \" set -eof pipefail; $UTILDIR/blast_outfmt6_replace_trans_id_w_gene_symbol.pl  ref_annot.cdsplus.fa ref_annot.cdsplus.allvsall.outfmt6 | gzip > ref_annot.cdsplus.allvsall.outfmt6.genesym.gz\" ";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdsplus.allvsall.outfmt6.genesym.gz.ok"));

    # index the blast hits:
    $cmd = "$UTILDIR/index_blast_pairs.pl $output_dir/blast_pairs.idx ref_annot.cdsplus.allvsall.outfmt6.genesym.gz";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/blast_pairs.idx.ok"));
    
    # remove blast pairs between genes that physically overlap on the genome
    $cmd = "$UTILDIR/index_blast_pairs.remove_overlapping_genes.pl $output_dir";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/blast_pairs.idx.ovrem.ok"));
    
    
    ##################################
    # blast across full cDNA sequences

    # extract the cDNA sequences
    $cmd = "$UTILDIR/gtf_file_to_feature_seqs.pl $gtf_file $genome_fa_file cDNA > ref_annot.cdna.fa";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdna.fa.ok"));
    
    $cmd = "makeblastdb -in ref_annot.cdna.fa -dbtype nucl";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdna.fa.blidx.ok"));
    
    $cmd = "blastn -query ref_annot.cdna.fa -db ref_annot.cdna.fa -max_target_seqs 10000 -outfmt 6 -evalue 1e-10 -num_threads $CPU -dust no  > ref_annot.cdna.allvsall.outfmt6";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdna.allvsall.outfmt6.ok"));
    
    $cmd = "$UTILDIR/isoform_blast_gene_chr_conversion.pl --blast_outfmt6 ref_annot.cdna.allvsall.outfmt6 --gtf $gtf_file > ref_annot.cdna.allvsall.outfmt6.toGenes";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdna.allvsall.outfmt6.toGenes.ok"));

    $cmd = "sort -k2,2 -k7,7 ref_annot.cdna.allvsall.outfmt6.toGenes > ref_annot.cdna.allvsall.outfmt6.toGenes.sorted";
    $pipeliner->add_commands(new Command($cmd, "ref_annot.cdna.allvsall.outfmt6.toGenes.sorted.ok"));

    $cmd = "$UTILDIR/build_chr_gene_alignment_index.pl --blast_genes ref_annot.cdna.allvsall.outfmt6.toGenes.sorted  --out_prefix $output_dir/trans.blast.align_coords";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/trans.blast.align_coords.ok"));
    
    
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
        $pipeliner->add_commands(new Command($cmd, "$output_dir/ref_genome.fa.gmap.ok"));
    }

        
    $pipeliner->run();

    exit(0);
}


####
sub check_for_required_tools {
    my (@tools) = @_;

    my $missing_flag = 0;

    foreach my $tool (@tools) {
        my $path = `which $tool`;
        if ($path =~ /\w/) {
            print STDERR "-found $tool at $path\n";
        }
        else {
            print STDERR "- *** ERROR, cannot locate required tool in PATH setting: $tool ***\n";
            $missing_flag = 1;
        }
    }

    if ($missing_flag) {
        die "Error, missing at least one required tool. See error messages and perform required software installations before running.";
    }

}
    


        
    
