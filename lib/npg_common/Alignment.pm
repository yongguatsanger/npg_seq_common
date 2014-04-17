#########
# Author:        gq1
# Created:       2009-04-17
#

package npg_common::Alignment;

use Carp;
use English qw(-no_match_vars);
use File::Temp qw(tempdir);
use Parallel::ForkManager;
use Moose;

with qw( npg_common::roles::software_location );

our $VERSION = '0';

has 'bwa_options' => (
    is            => 'rw',
    isa           => 'Maybe[Str]',
);

sub bwa_align_se {
  my ($self,$args_ref) = @_;

  my $query_fastq  = $args_ref->{fastq};
  my $bam_out  = $args_ref->{bam_out};
  my $ref_root  = $args_ref->{ref_root};

  my $hit_sequence_ids     = $args_ref->{hit_sequence_ids};

  my $work_dir = tempdir(CLEANUP => 1);
  my $sai_out = $work_dir.q{/sai};

  #alignment
  $self->bwa_align( $ref_root, $query_fastq, $sai_out );

  #sai to sam or bam
  $args_ref ->{sai} = $sai_out;

  if($hit_sequence_ids){
     #convert sai to bam and at the same time parsing each sam line
     $self->sai2bam_parsing_se($args_ref);
  }else{
     #convert sai to sam
     $self->sai2sam_se($args_ref);
  }

  #remove sai file when sam file generated to save disk space
  unlink $sai_out;

  return 1;
}

sub bwa_align_pe {
  my ( $self, $args_ref ) = @_;

  my $ref_root   = $args_ref->{ref_root};
  my $fastq1     = $args_ref->{fastq1};
  my $fastq2     = $args_ref->{fastq2};
  my $sam_out    = $args_ref->{sam_out};
  my $fork_align = $args_ref->{fork_align};

  my $work_dir = tempdir(CLEANUP => 1);

  my $sai_out_1 = $work_dir.q{/1.sai};
  my $sai_out_2 = $work_dir.q{/2.sai};

  #alignments
  $self->bwa_align_list( $ref_root, [$fastq1, $fastq2], [$sai_out_1, $sai_out_2], $fork_align);

  #sai to sam
  $args_ref->{sai_1} = $sai_out_1;
  $args_ref->{sai_2} = $sai_out_2;
  $self->sai2sam_pe($args_ref);

  #remove sai files
  unlink $sai_out_1;
  unlink $sai_out_2;

  return 1;
}

sub bwa_align {
  my ( $self, $ref_root, $query_fastq, $sai_out ) = @_;

  my $bwa_align_cmd = $self->_construct_bwa_aln_cmd( $ref_root, $query_fastq, $sai_out );

  print {*STDERR} ">Doing bwa alignment:\n $bwa_align_cmd\n" or croak $ERRNO;

  my $bwa_align_rs = system $bwa_align_cmd;
  if($bwa_align_rs != 0){
    croak "bwa align failed: $bwa_align_cmd";
  }

  return $bwa_align_cmd;
}

sub _construct_bwa_aln_cmd {
   my ($self, $ref_root, $query_fastq, $sai_out ) = @_;

   my $bwa_align_cmd = $self->bwa_cmd().q( aln );

   if( $self->bwa_options() ){
     $bwa_align_cmd .= $self->bwa_options();
   }

   $bwa_align_cmd .= qq( $ref_root $query_fastq > $sai_out);

  return $bwa_align_cmd;
}

sub bwa_align_list {
  my ( $self, $ref_root, $query_fastqs, $sai_out_files, $fork_align ) = @_;

  my $count = 0;

  if( $fork_align ){

    my $pm = Parallel::ForkManager->new( scalar @{$query_fastqs} );

    foreach my $fastq (@{$query_fastqs}){

      $count++;
      $pm->start and next;
      $self->bwa_align( $ref_root, $fastq, $sai_out_files->[$count-1] );
      $pm->finish;
    }
    $pm->wait_all_children;
  }else{

    foreach my $fastq (@{$query_fastqs}){

      $self->bwa_align( $ref_root, $fastq, $sai_out_files->[$count] );
      $count++;
    }
  }

  my $count2 = 0;
  my $bwa_aln_cmd = q{};
  foreach my $fastq (@{$query_fastqs}){
    $bwa_aln_cmd .= $self->_construct_bwa_aln_cmd( $ref_root, $fastq, $sai_out_files->[$count2]).q{;};
    $count2++;
  }

  return $bwa_aln_cmd;
}

sub sai2sam_se {
  my ( $self, $args ) = @_;

  my $ref_root  = $args->{ref_root};
  my $sai       = $args->{sai};
  my $fastq     = $args->{fastq};
  my $sam_out   = $args->{bam_out};

  my $bwa_sam_cmd = $self->bwa_cmd().qq( samse $ref_root $sai $fastq > $sam_out );

  print {*STDERR} ">bwa converting single end sai files to sam:\n $bwa_sam_cmd\n" or croak $ERRNO;

  my $bwa_sam_rs = system $bwa_sam_cmd;
  if($bwa_sam_rs != 0){
    croak "bwa align failed: $bwa_sam_cmd";
  }

  return 1;
}

sub sai2sam_pe {
  my ( $self, $args ) = @_;

  my $ref_root  = $args->{ref_root};
  my $sai_1     = $args->{sai_1};
  my $sai_2     = $args->{sai_2};
  my $fastq1    = $args->{fastq1};
  my $fastq2    = $args->{fastq2};
  my $sam_out   = $args->{sam_out};

  my $bwa_sam_cmd = $self->bwa_cmd().qq( sampe $ref_root $sai_1 $sai_2 $fastq1 $fastq2 > $sam_out );

  print {*STDERR} ">bwa converting paired end sai files to sam:\n $bwa_sam_cmd\n" or croak $ERRNO;

  my $bwa_sam_rs = system $bwa_sam_cmd;
  if($bwa_sam_rs != 0){
    croak "bwa align failed: $bwa_sam_cmd";
  }

  return 1;
}

1;

__END__

=head1 NAME

npg_common::Alignment

=head1 VERSION

=head1 SYNOPSIS

  my $oBWAAlignment = npg_common::Alignment->new(bwa_cmd => q{bwa});
  $oBWAAlignment->bwa_align_se($query_fastq, $sam, $reference_bwa);
  $oBWAAlignment->bwa_align_pe(
    {fastq1=>$query_fastq1, fastq2=>$query_fastq2, sam_out=>$sam, fork_align=>1, ref_root=>$reference_bwa});

=head1 DESCRIPTION

Class to do bwa alignment based on the fastq file

=head1 SUBROUTINES/METHODS

=head2 bwa_options

=head2 bwa_align - given one fastq file and reference location, doing bwa alignment and generate sai file

=head2 bwa_align_list - given a list of fastq files and output sai file list, doing bwa alignment

=head2 bwa_align_se - given one fastq file and reference location, doing bwa alignment and output in sam format

=head2 bwa_align_pe - given two fastq files and reference location, doing bwa alignment and generate one sam file. The two alignments could be done in parallel.

=head2 sai2sam_se - bwa convert sai file to sam file

=head2 sai2sam_pe - given two paired read sai file, convert them to sam file

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Carp

=item English

=item File::Temp qw(tempdir)

=item Parallel::ForkManager

=item Moose

=item npg_common::roles::software_location

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi, E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2009 GRL, by Guoying Qi

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.

=cut
