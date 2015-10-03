#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/lib");
use Pipeliner;


my $CPU = 4;

my $usage = <<__EOUSAGE__;

##################################################################################
#
#  Required:
#
#  --genome_fa  <string>           genome fasta file
#
#  --cdna_fa <string>              cdna fasta file
#                                  Note: header format must be:
#                                     transcript_id(tab)gene_id(tab)gene_symbol
#
#  --gtf <string>                  transcript structure annotation
#                                  Note: can restrict to coding genes
#
#  --blast_pairs <string>          transcript blastn results 
#                                  in BLAST+ '-outfmt 6'  format and gzipped!
#                                  Note: gene symbols must replace transcript IDs.
#
#  Optional:
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

&GetOptions ( 'h' => \$help_flag,

              # required
              'genome_fa=s' => \$genome_fa_file,
              'cdna_fa=s' => \$cdna_fa_file,
              'gtf=s' => \$gtf_file,
              'blast_pairs=s' => \$blast_pairs_file,

              # optional
              'output_dir=s' => \$output_dir,
              'CPU=i' => \$CPU,
    );


if ($help_flag) {
    die $usage;
}

unless ($genome_fa_file && $cdna_fa_file && $gtf_file && $blast_pairs_file) {
    die $usage;
}

unless ($blast_pairs_file =~ /\.gz$/) {
    die "Error, file $blast_pairs_file must be gzipped-compressed";
}


$genome_fa_file = &Pipeliner::ensure_full_path($genome_fa_file);
$cdna_fa_file = &Pipeliner::ensure_full_path($cdna_fa_file);
$gtf_file = &Pipeliner::ensure_full_path($gtf_file);
$blast_pairs_file = &Pipeliner::ensure_full_path($blast_pairs_file);


my $UTILDIR = $FindBin::Bin . "/util";

if ($output_dir) {
    unless (-d $output_dir) {
        mkdir $output_dir or die "Error, cannot mkdir $output_dir";
    }
    chdir $output_dir or die "Error, cannot cd to $output_dir";
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

    # build star index
    my $cmd = "ln -s $genome_fa_file genome.fa";
    $pipeliner->add_commands(new Command($cmd, "genome.fa.ok"));
    
    my $star_index = "genome.fa.star.idx";
    $cmd = "STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index "
            . " --twopassMode Basic "
            . " --genomeFastaFiles genome.fa "
            . " --limitGenomeGenerateRAM 40419136213 "
            . " --sjdbGTFfile $gtf_file "
            . " --sjdbOverhang 100 ";
    
    $pipeliner->add_commands(new Command($cmd, "$star_index.ok"));

    # build GMAP index
    $cmd = "gmap_build -D . -d genome.fa.gmap -T . -k 13 genome.fa";
    $pipeliner->add_commands(new Command($cmd, "genome.fa.gmap.ok"));


    ##########################
    # Prep the cdna fasta file
    
    $cmd = "ln -s $cdna_fa_file ref_cdna.fasta";
    $pipeliner->add_commands(new Command($cmd, "ref_cdna.fasta.ok"));

    # index the fasta file
    $cmd = "$UTILDIR/index_cdna_seqs.pl ref_cdna.fasta";
    $pipeliner->add_commands(new Command($cmd, "ref_cdna.fasta.idx.ok"));

    
    #######################
    # index the blast pairs
    
    $cmd = "$UTILDIR/index_blast_pairs.pl $blast_pairs_file blast_pairs.idx";
    $pipeliner->add_commands(new Command($cmd, "blast_pairs.idx.ok"));


    $pipeliner->run();

    exit(0);
}




        
    
