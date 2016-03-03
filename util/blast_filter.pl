#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 
use TiedHash;
use Data::Dumper;


my $Evalue = 1e-3;

my $usage = <<__EOUSAGE__;

###################################################################################################
#
# Required:
#
#  --fusion_preds <string>        preliminary fusion predictions
#                                 Required formatting is:  
#                                 geneA--geneB (tab) junction_read_count (tab) spanning_read_count (tab) ... rest
#
#  --genome_lib_dir <string>      genome lib directory
#
#
# Optional: 
##
#  -E <float>                     E-value threshold for blast searches (default: $Evalue)
#
####################################################################################################

This script will generate two files:

    \${fusion_preds}.post_blast_filter  : contains results filtered from likely false positive fusions w/ similar sequences.

    \${fusion_preds}.post_blast_filter.info : contains all input fusion predictions and those having blast-matches are annotated accordingly.


__EOUSAGE__

    ;

my $help_flag;

my $fusion_preds_file;
my $genome_lib_dir;

&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,
              
              'E=f' => \$Evalue,
                            
              'genome_lib_dir=s' => \$genome_lib_dir,
                            
    );


if (@ARGV) {
    die "Error, dont recognize arguments: @ARGV";
}


if ($help_flag) {
    die $usage;
}

unless ($fusion_preds_file && $genome_lib_dir) {
    die $usage;
}


my $BLAST_PAIRS_IDX;
my $blast_pairs_idx_file = "$genome_lib_dir/blast_pairs.idx";
if (-s $blast_pairs_idx_file) {
    $BLAST_PAIRS_IDX = new TiedHash( { use => $blast_pairs_idx_file } );
}
else {
    die "Error: cannot locate $blast_pairs_idx_file";
}

my %BLAST_CACHE;



=input_format:

0       ETV6--NTRK3
1       84
2       18
3       ONLY_REF_SPLICE
4       ETV6^ENSG00000139083.6
5       chr12:12022903:+
6       NTRK3^ENSG00000140538.12
7       chr15:88483984:-
8       comma-delim list of junction reads
9       comma-delim list of spanning frags

=cut


main: {

    my $final_preds_file = "$fusion_preds_file.post_blast_filter";
    open (my $final_ofh, ">$final_preds_file") or die "Error, cannot write to $final_preds_file";
    
    my $filter_info_file = "$fusion_preds_file.post_blast_filter.info";
    open (my $filter_ofh, ">$filter_info_file") or die "Error, cannot write to $filter_info_file";

    my %already_approved;

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    my $header = <$fh>;
    unless ($header =~ /^\#/) {
        die "Error, file $fusion_preds_file doesn't begin with a header line";
    }
    print $filter_ofh $header;
    print $final_ofh $header;
    
    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $fusion_name = $x[0];

        my ($geneA, $geneB) = split(/--/, $fusion_name);

        my @blast_info;

        unless ($already_approved{$geneA}->{$geneB}) {
            
            @blast_info = &examine_seq_similarity($geneA, $geneB);
            if (@blast_info) {
                push (@blast_info, "SEQ_SIMILAR_PAIR");
            }
        }
        
                
        if (@blast_info) {
            
            print $filter_ofh "#" . "$line\t" . join("\t", @blast_info) . "\n";
        }
        else {
            
            print $final_ofh "$line\n";
            print $filter_ofh "$line\n";
            
            $already_approved{$geneA}->{$geneB} = 1;
        }
        
    }
    close $fh;

    close $filter_ofh;
    close $final_ofh;
    
    

    exit(0);
}


####
sub examine_seq_similarity {
    my ($geneA, $geneB) = @_;
    
    #print STDERR "-examining seq similarity between $geneA and $geneB\n";

    my @blast_hits;

    # use pre-computed blast pair data
    if (my $hit = $BLAST_PAIRS_IDX->get_value("$geneA--$geneB")) {
        return($hit);
    }
    elsif ($hit = $BLAST_PAIRS_IDX->get_value("$geneB--$geneA")) {
        return($hit);
    }
    else {
        return();
    }
}
