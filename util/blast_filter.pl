#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 
use TiedHash;
use Data::Dumper;
use Gene_overlap_check;
use DelimParser;


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
####################################################################################################

This script will generate two files:

    \${fusion_preds}.post_blast_filter  : contains results filtered from likely false positive fusions w/ similar sequences.

    \${fusion_preds}.post_blast_filter.info : contains all input fusion predictions and those having blast-matches are annotated accordingly.


__EOUSAGE__

    ;

my $help_flag;

my $fusion_preds_file;
my $genome_lib_dir;

my $EXCLUDE_LOCI_OVERLAP_CHECK = 0;


&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,
              
              'genome_lib_dir=s' => \$genome_lib_dir,
                            
              'exclude_loci_overlap_check' => \$EXCLUDE_LOCI_OVERLAP_CHECK,
              
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



main: {

    my $final_preds_file = "$fusion_preds_file.post_blast_filter";
    open (my $final_ofh, ">$final_preds_file") or die "Error, cannot write to $final_preds_file";
    
    my $filter_info_file = "$fusion_preds_file.post_blast_filter.info";
    open (my $filter_ofh, ">$filter_info_file") or die "Error, cannot write to $filter_info_file";

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    my @column_headings = $delim_parser->get_column_headers();

    my $final_ofh_writer = new DelimParser::Writer($final_ofh, "\t", \@column_headings);
    my $filter_ofh_writer = new DelimParser::Writer($filter_ofh, "\t", [@column_headings, "FilterReason"]);
    
    
    
    my @fusions;
    
    while (my $row = $delim_parser->get_row()) { 
        my $fusion_name = $row->{'#FusionName'};
        
        my ($geneA, $geneB) = split(/--/, $fusion_name);

        my $J = 0;
        my $S = 0;
        if (defined $row->{JunctionReadCount}) {
        
            $J = ( (!defined($row->{est_J})) || $row->{est_J} eq "NA") ? $row->{JunctionReadCount} : $row->{est_J};
            $S = ( (!defined($row->{est_S})) ||  $row->{est_S} eq "NA") ? $row->{SpanningFragCount} : $row->{est_S};
        }
        
        my $num_LR = $row->{num_LR} || $row->{num_reads} || 0;
        
        if ( $J !~ /\d/) { $J = 0; }
        if ( $S !~ /\d/) { $S = 0; }
        
        if ( $num_LR eq "NA") { $num_LR = 0; }
        
        my $score = $num_LR + 4*$J + $S;
        
        push (@fusions, { row => $row,
                          fusion_name => $fusion_name,
                          geneA => $geneA,
                          geneB => $geneB,
                          score => $score,
              } );
        
    }
    close $fh;

    
    ## sort in order of score, descendingly
    @fusions = reverse sort {$a->{score}<=>$b->{score}
                             ||
                                 $b->{fusion_name} cmp $a->{fusion_name}  # for more stable sorting.
                             
                             

    } @fusions;
    

    ########################
    ##  Filter and report ##
    ########################
    
    #######################################################################
    ## Screen out paralog hits such that fusion A--B exists w/ a high score
    ##    and A--C exists with a lower score
    ##    and B has sequence similarity to C
    ##   in which case, we keep A--B and discard A--C
    ##   (unless B and C physically overlap on the genome!)
    #######################################################################

    my %AtoB;
    my %BtoA;
    
    my %already_approved;

    my $gene_overlap_checker;
    unless ($EXCLUDE_LOCI_OVERLAP_CHECK) {
        $gene_overlap_checker = new Gene_overlap_check("$genome_lib_dir/ref_annot.gtf.gene_spans");
    }
    
    foreach my $fusion (@fusions) {
        
        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};

	my $row = $fusion->{row};
        
        my @blast_info;

        if ($already_approved{$geneA}->{$geneB}) {
            
            # no op
        }
        else {
            
            @blast_info = &examine_seq_similarity($geneA, $geneB);
            if (@blast_info) {
    
                # immediately discarding seq-similar fusion partners
                push (@blast_info, "SEQ_SIMILAR_PAIR");
            }
            else {
                
                ## See if we already captured a fusion containing a paralog of the partner here:
                
                my $altB_href = $AtoB{$geneA};
                if ($altB_href) {
                    foreach my $altB (keys %$altB_href) {
                        my @blast = &examine_seq_similarity($geneB, $altB);
                        my $overlapping_genes_flag = ($EXCLUDE_LOCI_OVERLAP_CHECK) ? 0 : $gene_overlap_checker->are_genes_overlapping($geneB, $altB);
                        if (@blast && ! $overlapping_genes_flag) {
                            push (@blast, "ALREADY_EXAMINED:$geneA--$altB");
                            push (@blast_info, @blast);
                        }
                    }
                }
                
                my $altA_href = $BtoA{$geneB};
                if ($altA_href) {
                    foreach my $altA (keys %$altA_href) {
                        my @blast = &examine_seq_similarity($altA, $geneA);
                        my $overlapping_genes_flag = ($EXCLUDE_LOCI_OVERLAP_CHECK) ? 0 : $gene_overlap_checker->are_genes_overlapping($altA, $geneA);
                        if (@blast && ! $overlapping_genes_flag) {
                            push (@blast, "ALREADY_EXAMINED:$altA--$geneB");
                            push (@blast_info, @blast);
                        }
                    }
                }
            }
        }
        
        my $line = $fusion->{line};
        
        if (@blast_info) {

	    $row->{FilterReason} = join("^^^", @blast_info);
	    
	    $filter_ofh_writer->write_row($row);
	    
        }
        else {
            
            $final_ofh_writer->write_row($row);
            
            $already_approved{$geneA}->{$geneB} = 1;
        }

        $AtoB{$geneA}->{$geneB} = 1;
        $BtoA{$geneB}->{$geneA} = 1;
        
    }


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

