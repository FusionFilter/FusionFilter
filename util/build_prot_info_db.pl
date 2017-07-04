#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use JSON::XS;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Process_cmd;

my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --gtf <string>        : reference annotation for coding genes in gtf format
#                                      
#  --genome_fa <string>       : genome fasta file
#
#  --out_db_dir <string>      : directory to store prot_info.db.idx
#
##########################################################################


 
__EOUSAGE__

    ;


my $help_flag;
my $gtf_file;
my $genome_fa;
my $out_db_dir;

&GetOptions ( 'h' => \$help_flag,
              'gtf=s' => \$gtf_file,
              'genome_fa=s' => \$genome_fa,
              'out_db_dir=s' => \$out_db_dir,
    );

unless ($gtf_file && $genome_fa && $out_db_dir) {
    die $usage;
}


main: {

    unless (-d $out_db_dir) {
        &process_cmd("mkdir -p $out_db_dir");
    }
    my $prot_info_db_idx = "$out_db_dir/prot_info_db.idx";
    if (-e $prot_info_db_idx) {
        die "Error, $prot_info_db_idx already exists.  Please remove it or rename it before proceeding";
    }

    my $annot_manager = Annotation_manager->new($gtf_file, $genome_fa);

    $annot_manager->build_prot_info_db($prot_info_db_idx);
    
    
    exit(0);
}

        

#############################################################################
package Annotation_manager;
use strict;
use warnings;
use Gene_obj;
use Gene_obj_indexer;
use GTF_utils;
use Fasta_reader;
use TiedHash;
use Carp;

####
sub new {
    my ($packagename, $gtf_file, $genome_fa) = @_;

    my $self = {
        gene_to_CDS_features => {},  # gene_id -> cds_id -> cds_feature_obj
    };
    
    bless($self, $packagename);

    $self->parse_GTF_instantiate_featureset($gtf_file, $genome_fa);
    
    return($self);
}


sub TO_JSON {
    return { %{ shift() } };
}


####
sub get_gene_list {
    my ($self) = @_;

    my @gene_ids = keys %{$self->{gene_to_CDS_features}};

    return(@gene_ids);
}

####
sub get_CDS_features_for_gene {
    my $self = shift;
    my ($gene_id) = @_;

    my $cds_features_href = $self->{gene_to_CDS_features}->{$gene_id};
    
    if (ref $cds_features_href) {
        return(values %$cds_features_href);
    }

    else {
        return();
    }
}

####
sub toString {
    my ($self) = shift;
    
    my @gene_ids = $self->get_gene_list();

    my $text = "";
    
    foreach my $gene_id (@gene_ids) {
        my @cds_features = $self->get_CDS_features_for_gene($gene_id);
        foreach my $cds_feature (@cds_features) {
            $text .= $cds_feature->toString() . "\n";
        }
    }
    
    return ($text);
}
    
####
sub add_cds_and_pfam {
    my ($self) = shift;
    my ($cds_seqs_href, $pfam_hits_href) = @_;
    
    
    my @gene_ids = $self->get_gene_list();
    
    foreach my $gene_id (@gene_ids) {
        my @cds_features = $self->get_CDS_features_for_gene($gene_id);
        foreach my $cds_feature (@cds_features) {

            my $cds_id = $cds_feature->{cds_id};
            my $cds_seq = $cds_seqs_href->{$cds_id};
            unless ($cds_seq) {
                print STDERR "WARNING, no CDS sequence for $cds_id\n";
                $cds_seq = "";
            }
            
            $cds_feature->set_CDS_sequence($cds_seq);

            my $pfam_hits = $pfam_hits_href->{$cds_id};
            if (ref $pfam_hits) {
                $cds_feature->add_pfam_hits(@$pfam_hits);
            }
        }
    }

    return;
}

####
sub build_prot_info_db {
    my ($self) = shift;
    my ($prot_db_idx_file) = @_;
    
    my $tied_hash = new TiedHash( { create => $prot_db_idx_file } );
    
    
    my @gene_ids = $self->get_gene_list();
    
    my $coder = JSON::XS->new->convert_blessed;
    
    foreach my $gene_id (@gene_ids) {
        my @cds_features = $self->get_CDS_features_for_gene($gene_id);
        
        my $gene_id_store = $gene_id;
        $gene_id_store =~ s/\|.*$//;

        my $json = $coder->pretty->encode(\@cds_features);

        $tied_hash->store_key_value($gene_id_store, $json);
    }
    
    $tied_hash = undef; # closes it.
    
    return;
}


            
####
sub parse_GTF_instantiate_featureset {
    my ($self) = shift;
    my ($gtf_file, $genome_fa) = @_;

    my $gene_obj_indexer = {};

    ## associate gene identifiers with contig id's.
    &GTF_utils::index_GTF_gene_objs($gtf_file, $gene_obj_indexer);


    my $fasta_reader = new Fasta_reader($genome_fa);
    my %genome = $fasta_reader->retrieve_all_seqs_hash();

    my @all_gene_ids = keys %$gene_obj_indexer;
    my %contig_to_gene_list;
    foreach my $gene_id (@all_gene_ids) {
        my $gene_obj = $gene_obj_indexer->{$gene_id};

        my $contig = $gene_obj->{asmbl_id}
        or croak "Error, can't find contig id attached to gene_obj($gene_id) as asmbl_id val\n"
            . $gene_obj->toString();


        my $gene_list_aref = $contig_to_gene_list{$contig};
        unless (ref $gene_list_aref) {
            $gene_list_aref = $contig_to_gene_list{$contig} = [];
        }

        push (@$gene_list_aref, $gene_id);

    }

    foreach my $asmbl_id (sort keys %contig_to_gene_list) {

        my $genome_seq = $genome{$asmbl_id} or croak "Error, no sequence stored for contig: $asmbl_id";

        my @gene_ids = @{$contig_to_gene_list{$asmbl_id}};

        foreach my $gene_id (@gene_ids) {
            my $gene_obj_ref = $gene_obj_indexer->{$gene_id};
            
            $gene_obj_ref->create_all_sequence_types(\$genome_seq);
            
            foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {

                unless ($isoform->is_coding_gene()) { next; }

                $isoform->set_CDS_phases(\$genome_seq);
                
                my $isoform_id = $isoform->{Model_feat_name};
                my $gene_id = $isoform->{TU_feat_name};
                my $orient = $isoform->get_orientation();
                
                my $cds_seq = $isoform->get_CDS_sequence();

                my $cds_feature_obj = $self->get_CDS_feature($gene_id, $isoform_id);
                
                my @exons = $isoform->get_exons();
                foreach my $exon (@exons) {
                    if (my $cds_obj = $exon->get_CDS_obj()) {
                        
                        my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds_obj->get_coords();
                        my $phase = $cds_obj->get_phase();
                        
                        $cds_feature_obj->add_segment($asmbl_id, $cds_lend, $cds_rend, $orient, $phase);
                    }
                }
                
                $cds_feature_obj->set_CDS_sequence($cds_seq);
                $cds_feature_obj->refine();
            }
        }
    }
    
    return;
}


####
sub get_CDS_feature {
    my ($self) = shift;
    my ($gene_id, $cds_id) = @_;

    my $cds_feature_obj = $self->{gene_to_CDS_features}->{$gene_id}->{$cds_id};

    unless (ref $cds_feature_obj) {
        $cds_feature_obj = $self->{gene_to_CDS_features}->{$gene_id}->{$cds_id} = CDS_feature->new($gene_id, $cds_id);
    }

    return($cds_feature_obj);
}


####
sub _refine_CDS_features {
    my $self = shift;

    my @gene_ids = $self->get_gene_list();

    foreach my $gene_id (@gene_ids) {
        
        my $cds_features_href = $self->{gene_to_CDS_features}->{$gene_id};
        my @cds_feature_objs = values %$cds_features_href;

        foreach my $cds_feature_obj (@cds_feature_objs) {
            $cds_feature_obj->refine();
        }
    }
    
    return;
}



####################################################################################
package CDS_feature;
use strict;
use warnings;

####
sub new {
    my ($packagename) = shift;
    my ($gene_id, $cds_id) = @_;

    
    my $self = {
        phased_segments => [],
        
        gene_id => $gene_id,
        cds_id => $cds_id,
        
        refined_flag => 0,

        pfam_hits => [],

        cds_seq => "",
        
        
    };

    bless($self, $packagename);

    return($self);
}


sub TO_JSON {
    return { %{ shift() } };
}




####
sub add_segment {
    my ($self) = shift;
    my ($chr, $lend, $rend, $orient, $phase_beg) = @_;

    my $phased_segment = { chr => $chr,
                           
                           lend => $lend,
                           rend => $rend,
                           orient => $orient,
                           phase_beg => $phase_beg,

                           rel_lend => undef,
                           rel_rend => undef,
                           phase_end => undef,  ## all set on init
    };

    push (@{$self->{phased_segments}}, $phased_segment);
    
    return;
}

####
sub get_segments {
    my ($self) = shift;

    return(@{$self->{phased_segments}});

}


####
sub refine {
    my ($self) = shift;
    
    my @segments = $self->get_segments();

    @segments = sort {$a->{lend} <=> $b->{lend}} @segments;
    
    my $orient = $segments[0]->{orient};

    if ($orient eq '-') {
        @segments = reverse @segments;
    }

    my $sum_segs_len = 0;
    foreach my $segment (@segments) {

        my $seg_len = $segment->{rend} - $segment->{lend} + 1;
        my $phase_beg = $segment->{phase_beg};


        
        my $rel_lend = $sum_segs_len + 1;
        my $rel_rend = $sum_segs_len + $seg_len;

        my $phase_end = ".";
        if ($phase_beg ne ".") {
            my $adj_seg_len = $seg_len;
            $adj_seg_len += $phase_beg;
        
            $phase_end = ($adj_seg_len -1)  % 3;
        }
        
        $segment->{rel_lend} = $rel_lend;
        $segment->{rel_rend} = $rel_rend;
        $segment->{phase_end} = $phase_end;
        
        $sum_segs_len += $seg_len;
    
    }

    $self->{refined_flag} = 1;

    return;
}


####
sub toString {
    my ($self) = shift;
    
    my @segments = $self->get_segments();
    
    @segments = sort {$a->{lend} <=> $b->{lend}} @segments;

    my $orient = $segments[0]->{orient};
    if ($orient eq '-') {
        @segments = reverse @segments;
    }


    my $ret_text = "";
    
    foreach my $segment (@segments) {
        
        $ret_text .= join("\t", $self->{gene_id}, $self->{cds_id}, 
                          $segment->{lend}, $segment->{rend},
                          $segment->{orient}, 
                          $segment->{rel_lend}, $segment->{rel_rend},
                          $segment->{phase_beg}, $segment->{phase_end}) . "\n";
    }
    
    
    
    return ($ret_text);
        
}

####
sub set_CDS_sequence {
    my ($self) = shift;
    my ($cds_seq) = @_;

    $self->{cds_seq} = $cds_seq;

    return;
}

