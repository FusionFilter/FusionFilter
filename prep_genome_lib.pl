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

my $usage = <<__EOUSAGE__;

##################################################################################
#
#  Required by: STAR-Fusion, GMAP-Fusion, and FusionInspector
#
#  --genome_fa  <string>           genome fasta file
#
#  --gtf <string>                  transcript structure annotation
#                                     Note: can restrict to coding genes
#
#  --blast_pairs <string>          transcript blastn results 
#                                  in BLAST+ '-outfmt 6'  format and gzipped!
#                                     Note: gene symbols must replace transcript IDs.
#
# Required by: GMAP-Fusion and FusionInspector
#
#  --cdna_fa <string>              cdna fasta file
#                                     Note: header format must be:
#                                     transcript_id(tab)gene_id(tab)gene_symbol
#
# Required by: FusionInspector:
#
#  --count_kmers                   flag to include additional build steps required by count_kmers
#                                     (requires 'jellyfish' software exist in your PATH setting)
#
#  Misc options:
#
#  --output_dir <string>           output directory
#
#  --CPU <int>                     number of threads (defalt: $CPU)
#
##################################################################################



__EOUSAGE__

    ;


my $help_flag;
my $genome_fa_file;
my $cdna_fa_file;
my $gtf_file;
my $blast_pairs_file;
my $output_dir;
my $count_kmers;

&GetOptions ( 'h' => \$help_flag,

              # required for STAR-Fusion, FusionInspector, GMAP-fusion
              'genome_fa=s' => \$genome_fa_file,
              'gtf=s' => \$gtf_file,        
              
              # required for STAR-Fusion, FusionInspector, and GMAP-Fusion
              'blast_pairs=s' => \$blast_pairs_file,

              # required for GMAP-fusion and FusionInspector
              'cdna_fa=s' => \$cdna_fa_file,
              
              # optional
              'output_dir=s' => \$output_dir,
              'CPU=i' => \$CPU,
    
              # required for FusionInspector
              'count_kmers' => \$count_kmersPrep,
    );


if ($help_flag) {
    die $usage;
}

unless ($genome_fa_file && $gtf_file && $blast_pairs_file) {
    die $usage;
}

unless ($blast_pairs_file =~ /\.gz$/) {
    die "Error, file $blast_pairs_file must be gzipped-compressed";
}


$genome_fa_file = &Pipeliner::ensure_full_path($genome_fa_file);
$cdna_fa_file = &Pipeliner::ensure_full_path($cdna_fa_file);
$gtf_file = &Pipeliner::ensure_full_path($gtf_file);
$blast_pairs_file = &Pipeliner::ensure_full_path($blast_pairs_file);
$output_dir = &Pipeliner::ensure_full_path($output_dir);

my $UTILDIR = $FindBin::Bin . "/util";

unless ($output_dir) {
    $output_dir = cwd();
}

my @tools_required = qw(STAR gmap_build bowtie-build);
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
    
    my $cmd;
    
    unless (-e "$output_dir/ref_genome.fa") {
        $cmd = "ln -s $genome_fa_file $output_dir/ref_genome.fa";
        $pipeliner->add_commands(new Command($cmd, "$output_dir/_ref_genome.fa.ok"));
    }

    # build star index
    my $star_index = "$output_dir/ref_genome.fa.star.idx";
    unless (-d $star_index) {
        mkpath $star_index or die "Error, cannot mkdir $star_index";
    }

    
    $cmd = "STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index "
            . " --twopassMode Basic "
            . " --genomeFastaFiles $output_dir/ref_genome.fa "
            . " --limitGenomeGenerateRAM 40419136213 "
            . " --sjdbGTFfile $gtf_file "
            . " --sjdbOverhang 100 ";
    
    $pipeliner->add_commands(new Command($cmd, "$star_index.ok"));

    # build GMAP index
    $cmd = "gmap_build -D $output_dir -d ref_genome.fa.gmap -k 13 $output_dir/ref_genome.fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir/_ref_genome.fa.gmap.ok"));


    ###############################
    ## symlink the annotation file
    
    unless (-e "$output_dir/ref_annot.gtf") {
        $cmd = "ln -sf $gtf_file $output_dir/ref_annot.gtf";
        $pipeliner->add_commands(new Command($cmd, "ref_annot.gtf.ok"));
    }
    

    #######################
    # index the blast pairs
    
    $cmd = "$UTILDIR/index_blast_pairs.pl $blast_pairs_file $output_dir/blast_pairs.idx";
    $pipeliner->add_commands(new Command($cmd, "$output_dir/_blast_pairs.idx.ok"));
    

    ##########################
    # Prep the cdna fasta file
    
    if ($cdna_fa_file) {
        unless (-e "$output_dir/ref_cdna.fasta") {
            $cmd = "ln -s $cdna_fa_file $output_dir/ref_cdna.fasta";
            $pipeliner->add_commands(new Command($cmd, "$output_dir/_ref_cdna.fasta.ok"));
        }
                
        # index the fasta file
        $cmd = "$UTILDIR/index_cdna_seqs.pl $output_dir/ref_cdna.fasta";
        $pipeliner->add_commands(new Command($cmd, "$output_dir/_ref_cdna.fasta.idx.ok"));
        
        # build the bowtie index:
        $cmd = "bowtie-build $output_dir/ref_cdna.fasta $output_dir/ref_cdna.fasta";
        $pipeliner->add_commands(new Command($cmd, "$output_dir/_ref_cdna.fasta.bowtie_idx.ok"));
    }
            
    if ($count_kmers) {
        my $cmd = "jellyfish count -t $CPU -m 25 -s 1000000000 --canonical $output_dir/ref_genome.fa";
        $pipeliner->add_commands(new Command($cmd, "$output_dir/_jf.count.ok"));

        $cmd = "jellyfish dump -L 2 mer_counts.jf > $output_dir/ref_genome.jf.min2.kmers";
        $pipeliner->add_commands(new Command($cmd, "$output_dir/_jf.dump.ok"));

        $cmd = "rm mer_counts.jf";
        $pipeliner->add_commands(new Command($cmd, "$output_dir/_jf.cleanup.ok"));
        
    }
    
    $pipeliner->run();

    exit(0);
}




        
    
