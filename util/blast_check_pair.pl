#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 
use TiedHash;
use Data::Dumper;


my $genome_lib_dir = $ENV{CTAT_GENOME_LIB} || die "Error, need CTAT_GENOME_LIB set";


my $usage = "\n\n\tusage: $0 geneA geneB\n\n";

my $geneA = $ARGV[0] or die $usage;
my $geneB = $ARGV[1] or die $usage;

my $BLAST_PAIRS_IDX;
my $blast_pairs_idx_file = "$genome_lib_dir/blast_pairs.idx";
if (-s $blast_pairs_idx_file) {
    $BLAST_PAIRS_IDX = new TiedHash( { use => $blast_pairs_idx_file } );
}
else {
    die "Error: cannot locate $blast_pairs_idx_file";
}



my @blast_info = &examine_seq_similarity($geneA, $geneB);
if (@blast_info) {
    print "@blast_info\n";
}
else {
    print "no hits.\n";
    
}

exit(0);


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

