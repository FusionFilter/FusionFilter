#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use TiedHash;
use Carp;
use PerlIO::gzip;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE;

###########################################################
#
# --blast_outfmt6 <string>      blast pairs 
#
# --gtf <string>                genome annotation gtf file
#
# --out_prefix                  output file prefix
#
###########################################################

__EOUSAGE

    ;

my $blast_file;
my $gtf_file;
my $out_prefix;


my $help_flag;


&GetOptions ( 'h' => \$help_flag,
              'blast_outfmt6=s' => \$blast_file,
              'gtf_file=s' => \$gtf_file,
              'out_prefix=s' => \$out_prefix,
    );


unless ($blast_file && $gtf_file && $out_prefix) {
    die $usage;
}


main: {

    my %isoform_info = &parse_isoform_gtf($gtf_file);

    #my $idx = new TiedHash( { create => $output_index_filename } );

    

    

    exit(0);
}


####
sub parse_isoform_gtf {
    my ($gtf_file) = @_;

    my %isoform_info;
    
    open(my $fh, $gtf_file) or die "Error, cannot open file: $gtf_file";
    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $chr = $x[0];
        my $feat_type = $x[2];

        unless($feat_type eq 'exon') { next; }
        
        
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        
        my $info = $x[8];

        my $gene_name;
        if ($info =~ /gene_name \"([^\"]+)/) {
            $gene_name = $1;
        }
        elsif ($info =~ /gene_id \"([^\"]+)/) {
            $gene_name = $1;
        }
        else {
            print STDERR "-not finding gene_id or gene_name for entry: $line\nskipping...\n";
            next;
        }

        my $transcript_id;
        if ($info =~ /transcript_id \"([^\"]+)/) {
            $transcript_id = $1;
        }
        else {
            print STDERR "-not finding transcript_id for line $line\nskipping...\n";
            next;
        }
        

        
        my $struct = $isoform_info{$transcript_id};
        unless ($struct) {
            $struct = $isoform_info{$transcript_id} = { gene_id => $gene_name,
                                                        transcript_id => $transcript_id,
                                                        chr => $chr,
                                                        orient => $orient,
                                                        exons => [],
            };
        }

        push (@{$struct->{exons}}, { lend => $lend,
                                     rend => $rend,
              } );

    }
    close $fh;
    

    foreach my $struct (values (%isoform_info)) {
        &set_rel_coords($struct);
    }

    return(%isoform_info);
}

####
sub set_rel_coords {
    my ($struct) = @_;

    my @exons = sort {$a->{lend} <=> $rend} @{$struct->{exons}};

    if ($struct->{orient} eq '-') {
        @exons = reverse @exons;
    }

    my $last_rel_pos = 0;
    
    foreach my $exon (@exons) {

        my ($lend, $rend) = ($exon->{lend}, $exon->{rend});
        
        my $exon_len = $rend - $lend + 1;

        $exon->{rel_lend} = $last_rel_pos + 1;
        $exon->{rel_rend} = $last_rel_pos + $exon_len;

        $last_rel_pos += $exon_len;
    }

    return;
}
    
