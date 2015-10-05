#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/lib");
use Fasta_reader;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 
use TiedHash;


my $Evalue = 1e-3;
my $tmpdir = "/tmp";

my $MAX_PROMISCUITY = 3;  # perhaps a poor choice of words, but still a best fit IMHO.


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
#  --tmpdir <string>              file for temporary files (default: $tmpdir)
#
#  --max_promiscuity <int>               maximum number of partners allowed for a given fusion. Default: $MAX_PROMISCUITY
#
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
              'tmpdir=s' => \$tmpdir,
              
              'max_promiscuity=i' => \$MAX_PROMISCUITY,
                   
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

    unless (-d $tmpdir) {
        mkdir $tmpdir or die "Error, cannot mkdir $tmpdir";
    }
    
    
    my $final_preds_file = "$out_prefix.final";
    open (my $final_ofh, ">$final_preds_file") or die "Error, cannot write to $final_preds_file";
    
    my $filter_info_file = "$fusion_preds_file.post_blast_n_promisc_filter";
    open (my $filter_ofh, ">$filter_info_file") or die "Error, cannot write to $filter_info_file";
    
    my @fusions;
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    my $header = <$fh>;
    unless ($header =~ /^\#/) {
        die "Error, file $fusion_preds_file doesn't begin with a header line";
    }
    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $fusion_name = $x[0];
        my $J = $x[1];
        my $S = $x[2];


        my ($geneA, $geneB) = split(/--/, $fusion_name);

        #my $score = sqrt($J**2 + $S**2);
        my $score = $J*4 + $S;
        
        
        my $fusion = { fusion_name => $fusion_name,
                       
                       J => $J,
                       S => $S,
                       
                       geneA => $geneA,
                       geneB => $geneB,
                       
                       score => $score, 
                                              
                       line => $line,
        };
    
        push (@fusions, $fusion); 
    }
    
    print $filter_ofh $header;
    print $final_ofh $header;
    
    @fusions = reverse sort {$a->{score} <=> $b->{score} } @fusions;
    
    @fusions = &remove_promiscuous_fusions(\@fusions, $filter_ofh, $MAX_PROMISCUITY);
        
    ########################
    ##  Filter and report ##
    ########################
    

    my %AtoB;
    my %BtoA;
    
    my %already_approved;

    foreach my $fusion (@fusions) {

        #my ($geneA, $geneB) = ($fusion->{geneA}, $fusion->{geneB});
        my ($geneA, $geneB) = split(/--/, $fusion->{fusion_name});


        #print STDERR "-encountered $geneA--$geneB\n\n";
        
        #print STDERR "Examining fusion: $geneA--$geneB\n";
        my @blast_info;

        if ($already_approved{$geneA}->{$geneB}) {
           
            #print STDERR "-already approved $geneA--$geneB\n";
        }
        else {
            
            @blast_info = &examine_seq_similarity($geneA, $geneB);
            if (@blast_info) {
                push (@blast_info, "SEQ_SIMILAR_PAIR");
            }
            else {
                
                my $altB_href = $AtoB{$geneA};
                if ($altB_href) {
                    foreach my $altB (keys %$altB_href) {
                        my @blast = &examine_seq_similarity($geneB, $altB);
                        if (@blast) {
                            push (@blast, "ALREADY_EXAMINED:$geneA--$altB");
                            push (@blast_info, @blast);
                        }
                    }
                }
                my $altA_href = $BtoA{$geneB};
                if ($altA_href) {
                    foreach my $altA (keys %$altA_href) {
                        my @blast = &examine_seq_similarity($altA, $geneA);
                        if (@blast) {
                            push (@blast, "ALREADY_EXAMINED:$altA--$geneB");
                            push (@blast_info, @blast);
                        }
                    }
                }
            }
        }
        
        my $line = $fusion->{line};
        
        if (@blast_info) {
            
            print $filter_ofh "#" . "$line\t" . join("\t", @blast_info) . "\n";
        }
        else {
            
            print $final_ofh "$line\n";
            print $filter_ofh "$line\n";
            
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
    

####
sub remove_promiscuous_fusions {
    my ($fusions_aref, $filter_ofh, $max_promiscuity) = @_;

    
    my %partners;

    foreach my $fusion (@$fusions_aref) {

        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};

        $partners{$geneA}->{$geneB}++;
        $partners{$geneB}->{$geneA}++;

    }

    my @ret_fusions;
    
    foreach my $fusion (@$fusions_aref) {

        my $geneA = $fusion->{geneA};
        my $geneB = $fusion->{geneB};

        my $num_geneA_partners = scalar(keys %{$partners{$geneA}});
        my $num_geneB_partners = scalar(keys %{$partners{$geneB}});

        if ($num_geneA_partners > $max_promiscuity || $num_geneB_partners > $max_promiscuity) {

            print $filter_ofh "#" . $fusion->{line} . "\tFILTERED DUE TO reached max promiscuity ($max_promiscuity), num_partners($geneA)=$num_geneA_partners and num_partners($geneB)=$num_geneB_partners\n";
        }
        else {
            push (@ret_fusions, $fusion);
        }
    }


    return(@ret_fusions);
}

