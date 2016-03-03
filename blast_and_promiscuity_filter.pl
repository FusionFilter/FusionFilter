#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 

my $Evalue = 1e-3;

my $MAX_PROMISCUITY = 3;  # perhaps a poor choice of words, but still a best fit IMHO.
my $MIN_PCT_DOM_PROM = 20;

my $usage = <<__EOUSAGE__;


###################################################################################################
#
# Required:
#
#  --fusion_preds <string>        preliminary fusion predictions
#                                 Required formatting is:  
#                                 geneA--geneB (tab) junction_read_count (tab) spanning_read_count (tab) ... rest
#
#  --out_prefix <string>          prefix for output filename (will tack on .final and .final.abridged)
#
#  --genome_lib_dir <string>      genome lib directory
#
# Optional: 
##
#  -E <float>                     E-value threshold for blast searches (default: $Evalue)
#
#  --max_promiscuity <int>               maximum number of partners allowed for a given fusion. Default: $MAX_PROMISCUITY
#
#  --min_pct_dom_promiscuity <int>       for promiscuous fusions, those with less than this support of the dominant scoring pair 
#                                        are filtered prior to applying the max_promiscuity filter.
#                                        (default: $MIN_PCT_DOM_PROM)
#
####################################################################################################


__EOUSAGE__

    ;

my $help_flag;

my $fusion_preds_file;
my $out_prefix;
my $genome_lib_dir;

&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,
              
              'out_prefix=s' => \$out_prefix,

              'E=f' => \$Evalue,
              
              'max_promiscuity=i' => \$MAX_PROMISCUITY,

              'min_pct_dom_promiscuity=i' => \$MIN_PCT_DOM_PROM,
                   
              'genome_lib_dir=s' => \$genome_lib_dir,
              
              
    );

if (@ARGV) {
    die "Error, dont recognize arguments: @ARGV";
}


if ($help_flag) {
    die $usage;
}

unless ($fusion_preds_file && $genome_lib_dir && $out_prefix) {
    die $usage;
}



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


my $UTILDIR = "$FindBin::Bin/util";


main: {

    ## run blast filter
    
    my $cmd = "$UTILDIR/blast_filter.pl --fusion_preds $fusion_preds_file -E $Evalue --genome_lib_dir $genome_lib_dir";
    &process_cmd($cmd);

    my $blast_filtered_preds = "$fusion_preds_file.post_blast_filter";
    unless (-s $blast_filtered_preds) {
       confess "Error, no output file: $blast_filtered_preds  was generated.";
    }

    ## run the promiscuity filter
    
    $cmd = "$UTILDIR/promiscuity_filter.pl --fusion_preds $blast_filtered_preds "
        . " --max_promiscuity $MAX_PROMISCUITY "
        . " --min_pct_dom_promiscuity $MIN_PCT_DOM_PROM ";

    &process_cmd($cmd);

    my $post_promisc_outfile = "$blast_filtered_preds.post_promisc_filter";
    unless (-s $post_promisc_outfile) {
        confess "Error, no output file: $post_promisc_outfile";
    }

    ## just copy it to the expected output file name:

    $cmd = "cp $post_promisc_outfile $out_prefix.final";
    &process_cmd($cmd);
    

    exit(0);
}



####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
        
    my $ret = system($cmd);
    if ($ret) {

        die "Error, cmd $cmd died with ret $ret";
    }
    
    return;
}
    
