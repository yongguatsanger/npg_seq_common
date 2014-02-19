#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2010 01 15
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::run::Bam;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use File::Basename;
use File::Spec;
use List::Util qw(first);
use IPC::Open3;
use JSON;
use Perl6::Slurp;

use st::api::lims;
use npg_common::irods::Loader;
use npg_common::sequence::BAM_MarkDuplicate;

with qw/
        MooseX::Getopt
        npg_tracking::illumina::run::short_info
        npg_tracking::illumina::run::folder
        npg_common::roles::log
        npg_common::roles::software_location
       /;
with qw{npg_tracking::illumina::run::long_info};
with qw{npg_common::roles::run::lane::tag_info};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $DEFAULT_ROOT_DIR       => q{/seq/};
Readonly::Scalar our $DEFAULT_LANE_NUMBERS   => 8;
Readonly::Scalar our $EXIT_CODE_SHIFT        => 8;


Readonly::Array  our @LIMS_PROPERTIES => qw/
                                            sample_name
                                            sample_public_name
                                            sample_accession_number
                                            sample_common_name
                                            sample_id
                                            study_id
                                            study_name
                                            study_title
                                            study_accession_number
                                          /;


## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::run::Bam

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head2 positions

a list of postions of the run

=cut
has 'positions'    => (isa           => q{ArrayRef},
                       is            => q{rw},
                       lazy_build    => 1,
                       documentation => q{a list of postions of the run},
                      );

sub _build_positions {
    my $self = shift;

    my $lims = st::api::lims->new(id_run => $self->id_run());
    my $positions = [ map {$_->position()} $lims->associated_child_lims ];

    return $positions;
}

=head2 alignment_filter_in_plex

to indicate spliked phix and non-consented reads split in plex level, not at lane level

=cut
has 'alignment_filter_in_plex'  => (isa           => q{Bool},
                                    is            => q{rw},
                                    lazy_build    => 1,
                                    documentation => q{indicate spliked phix and non-consented reads split in plex level, not at lane level},
                                   );

sub _build_alignment_filter_in_plex {
    my $self = shift;
    my @plex_filtered_files    = grep { m/ [#] \d+ _ (?:human|phix) [.] bam $/mxs } @{$self->file_list()};
    return scalar @plex_filtered_files  > 0;
}

=head2 nonhuman_not_in_filename

to indicate no nonhuman in target bam filename

=cut
has 'nonhuman_not_in_filename'  => (isa           => q{Bool},
                                    is            => q{rw},
                                    default       => 0,
                                    documentation => q{indicate no nonhuman in target bam filename},
                                   );

sub _build_run_folder {
    my ($self) = @_;
    if (! ($self->_given_path or $self->has_id_run or $self->has_name)){
       croak 'need a path to work out a run_folder';
    }
    return first {$_ ne q()} reverse File::Spec->splitdir($self->runfolder_path);
}

=head2 file_list

bam file list

=cut
has 'file_list'    => (isa           => q{ArrayRef},
                       is            => q{rw},
                       lazy_build    => 1,
                       documentation => q{bam file name list},
                      );
sub _build_file_list {
  my $self = shift;

  my @file_list = ();

  my $positions = $self->positions();
  my $positions_string = q{*};
  if( $positions && scalar @{$positions} > 0 ) {
    $positions_string = join q{,}, @{$positions};
    $positions_string = qq{{$positions_string}};
  }

  @file_list = glob  $self->archive_path().q{/}.$self->id_run.q{_}.$positions_string.q{*.bam};

  my @file_list_plex = glob $self->archive_path().q{/lane}.$positions_string.q{/*.bam};

  push @file_list, @file_list_plex;

  return \@file_list;
}

=head2 cal_file_list

cal file list

=cut
has 'cal_file_list'    => (isa           => q{ArrayRef},
                       is            => q{rw},
                       lazy_build    => 1,
                       documentation => q{cal file name list},
                      );
sub _build_cal_file_list {
  my $self = shift;

  my @file_list = ();

  @file_list = glob  $self->recalibrated_path().q{/}.q{*_purity_*};

  return \@file_list;
}

=head2 number_of_reads_list

number of reads for each bam file

=cut

has 'number_of_reads_list'    => (isa           => q{HashRef},
                                  is            => q{rw},
                                  default       => sub { {} },
                                  documentation => q{number of reads for each bam file},
                                );

=head2 number_of_reads_list_irods

number of reads for each bam file from irods meta data

=cut

has 'number_of_reads_list_irods' => (isa           => q{HashRef},
                                     is            => q{rw},
                                     default       => sub { {} },
                                     documentation => q{number of reads for each bam file from irods meta data},
                                   );

=head2 exclude_bam_index

a flag to include bam index file in the loading,
The index bam file must be in the same directory as bam file,
and the file name ends with .bai

=cut

has 'exclude_bam_index'   => (isa           => q{Bool},
                              is            => q{rw},
                              documentation => q{flag to exclue bam index file loading},
                             );

=head2 _bam_indexed_cache

a hash to cache bam indexed or not

=cut

has '_bam_indexed_cache'  => (isa           => q{HashRef},
                              is            => q{ro},
                              default       => sub { {} },
                             );
=head2 _bam_indexed_cache_irods

a hash to cache bam indexed or not based on irods meta data alignment

=cut

has '_bam_indexed_cache_irods'  => (isa           => q{HashRef},
                                    is            => q{ro},
                                    default       => sub { {} },
                                   );

=head2 alt_process
option non-standard process (eg 'casava')
=cut
has 'alt_process' => (isa          => 'Str',
                      is           => 'rw',
                      default      => q{},
                      documentation => 'non-standard process used',
                     );

=head2 collection
sub directory within irods to store results
=cut
has 'collection' => (isa            => 'Str',
                       is             => 'rw',
                       lazy_build     => 1,
                     documentation => 'collection within irods to store results',
                      );

sub _build_collection {
  my $self = shift;
  my $collection = $DEFAULT_ROOT_DIR . $self->id_run();
  if ($self->alt_process()) { $collection .= q{/} . $self->alt_process(); }
  return $collection;
}

=head2 process

main method to call, get each bam file, load them and add meta data

=cut

sub process{ ## no critic (Subroutines::ProhibitExcessComplexit)
  my $self = shift;

  my $id_run = $self->id_run();

  $self->log("Loading bam files for run $id_run");
  if( $self->positions() && ( scalar @{$self->positions()} ) ){
    $self->log( 'Lane '. (join q{, }, @{$self->positions}) );
  }


  # Now we handle the purity files
  my $cal_file_list = $self->cal_file_list();
  $self->log('There are '.( scalar @{ $cal_file_list} ). ' purity files to load');
  foreach my $file (@{ $cal_file_list}){

    $self->log("---- Loading file $file");

    my $file_name = basename($file);
    my ($lane)        = $file_name =~ m/ \d+ _ (\d) /mxs;

    my $meta  = {
                  id_run => $id_run,
                  lane   => $lane,
                  type   => q{purity},
                  alt_process => $self->alt_process(),
                };
    my $loader = npg_common::irods::Loader->new({
       file        => $file,
       collection  => $self->collection(),
       meta_data   => $meta,
    });

    $loader->run();
  }

  # Handle the BAM files
  my $file_list = $self->file_list();
  $self->log('There are '.( scalar @{ $file_list} ). ' bam files to load');
  my $file_count = 0;
  foreach my $file (@{ $file_list}){

    $self->log("---- Loading file $file");

    my $file_name = basename($file);
    my ($lane)        = $file_name =~ m/ \d+ _ (\d) /mxs;
    my ($tag_index)   = $file_name =~ m/ [#] (\d+) /mxs;
    my ($split) = $file_name =~ m/ _ (xahuman|yhuman|human|phix) /mxs;

    my $names = $self->get_study_library_sample_names($id_run, $lane, $tag_index);

    my $is_phix_control = $names->{is_phix_control};

    if( $is_phix_control && (! defined $tag_index) ) {
       $self->log('Phix control lane, IGNORED');
       next;
    }

    my $meta  = {
                  id_run => $id_run,
                  lane   => $lane,
                  type   => q{bam},
                  target => 1,
                };

    if(defined $tag_index){
       $meta->{tag_index} = $tag_index;

       my $tag_list = $self->tag_list($lane);
       my $tag = $tag_list->{$tag_index};
       if($tag){
         $meta->{tag} = $tag;
       }
       if($tag_index == 0){
         $meta->{target} = 0;
       }
    }

    if($split){
       $meta->{alignment_filter} = $split;
       if ($split ne 'yhuman') { $meta->{target} = 0 };
    }

    $meta->{is_paired_read} = $self->is_paired_read();

    my @irods_lims_meta_list = (@LIMS_PROPERTIES, q{sample}, q{study}, q{library}, q{library_id});

    foreach my $meta_name (@irods_lims_meta_list){
        my $values = $names->{$meta_name};
        if($values && scalar @{$values}){
           $meta->{$meta_name} = $values;
        }
    }

    my $reference_hashref = $self->get_bam_reference($file);
    my $reference = $reference_hashref->{reference};
    my $alignment = $reference_hashref->{alignment};
    if( $reference &&  $alignment){
       $meta->{reference} = $reference;
    }
    $meta->{alignment} = $alignment;

    my $total_reads = $self->get_number_of_reads($file);
    if( defined $total_reads ) {
       $meta->{total_reads} = $total_reads;
    }
    if( defined $total_reads && $total_reads == 0) {
       $meta->{alignment} = 0;
    }

    if ($self->alt_process()) {
        $meta->{alt_process} = $self->alt_process();
        if ($meta->{target}) {
            $meta->{target} = 0;
            $meta->{alt_target} = 1;
        }
    }

    my $loader = npg_common::irods::Loader->new({
       file        => $file,
       collection  => $self->collection(),
       meta_data   => $meta,
    });

    #give null permission for public group for a bam file with nonconsented sequence
    my @permissions = $split && ($split =~ /human/smx) && ($split ne 'yhuman') ? (q{null public}) : ();
    #if not nonconsented try to set read to group for study and then remove public read
    if( not @permissions  and my @study_id = @{$meta->{study_id}||[]}) {
      push @permissions, (map {qq{read ss_$_}} @study_id), q{null public};
    }
    $loader->chmod_permissions( \@permissions );

    if($loader->run()){
       $file_count++;
    }

    if( defined $total_reads && $total_reads != 0 && $alignment ) {
       $self->_load_bam_index_file($file, $id_run, \@permissions);
    }else {
       $self->log('No reads or no alignment in bam file so bam index file not loaded');
    }

	$self->_load_bam_other_file($file, $id_run, 0, '.bamcheck');
	$self->_load_bam_other_file($file, $id_run, 0, '.flagstat');
	$self->_load_bam_other_file($file, $id_run, \@permissions, '.cram');
	$self->_load_bam_other_file($file, $id_run, 0, '_quality_cycle_caltable.txt');
	$self->_load_bam_other_file($file, $id_run, 0, '_quality_cycle_surv.txt');
	$self->_load_bam_other_file($file, $id_run, 0, '_quality_error.txt');

  }

  $self->log("---- $file_count bam files loaded for run $id_run");

  if(!$self->bam_fully_archived){
     croak 'BAM files not fully archived';
  }else{
     $self->log('BAM files are fully archived');
  }

  return 1;
}

sub _load_bam_index_file {
    my ($self, $bam_file, $id_run, $permissions) = @_;

    if( $self->exclude_bam_index() ){
       return;
    }

    my $index_file = $bam_file;
    $index_file =~ s/bam$/bai/mxs;

    if( !-e $index_file ){
       if( !$self->_bam_no_alignment($bam_file) ){
          croak "Bam with alignment but index file does not exist: $index_file";
       }else{
          $self->log("Bam without alignment and index file does not exist: $index_file");
          return;
       }
    }

    $self->log("-Loading bam index file: $index_file");
    my $loader = npg_common::irods::Loader->new({
       file        => $index_file,
       collection  => $self->collection(),
       meta_data   => {type => q{bai}, alt_process => $self->alt_process()},
    });
    if( $permissions and @{$permissions} ){
       $loader->chmod_permissions($permissions);
    }
    $loader->run();

    return;
}

sub _load_bam_other_file {
	my ($self, $bam_file, $id_run, $permissions, $file_suffix) = @_;

	my $other_file = $bam_file;
	$other_file =~ s/[.]bam$/$file_suffix/mxs;

	if (! -e $other_file) {
		$self->log("Can't find file $other_file to load to iRods");
		return;
	}
	$self->log("-Loading other file: $other_file");
	my $loader = npg_common::irods::Loader->new({
		file=> $other_file,
		collection=> $self->collection(),
		meta_data=> {type => $file_suffix, alt_process => $self->alt_process(),},
	});
    if( $permissions and @{$permissions} ){
       $loader->chmod_permissions($permissions);
    }
    $loader->run();

    return;
}

=head2 get_study_library_sample_names

method to get study, library and sample names for a bam file

=cut

sub get_study_library_sample_names {
   my ($self, $id_run, $position, $tag_index) = @_;

   my $lims = st::api::lims->new(id_run => $id_run, position => $position, tag_index => $tag_index);

   my $spiked_phix_tag_index = $lims->spiked_phix_tag_index();
   my @all_lims = $lims->children; # in case we deal with tag 0
   if (!@all_lims) {
       @all_lims = ($lims);
   }

   my $contains_human   = $lims->contains_nonconsented_human;
   my $contains_xahuman = $lims->contains_nonconsented_xahuman;
   my $contains_yhuman = $lims->separate_y_chromosome_data;

   my $data = {};

   foreach my $prop (@LIMS_PROPERTIES) {
       $data->{$prop} = {};
   }
   foreach my $l (@all_lims) {

       my $l_tag_index = $l->tag_index();

       if( $spiked_phix_tag_index && (defined $l_tag_index) && $l_tag_index == $spiked_phix_tag_index && scalar @all_lims != 1){
          next;
       }
       foreach my $prop (@LIMS_PROPERTIES) {
	      if ($l->$prop) {
               $data->{$prop}->{$l->$prop} = 1;
	      }
       }
   }

   my $names = { study           => [sort keys %{$data->{study_name}}],
                 library         => [$lims->library_name || $lims->library_id],
                 library_id      => [$lims->library_id],
                 sample          => [sort keys %{$data->{sample_name}}],
                 is_phix_control => $lims->is_control,
                 is_human_split  => $contains_human,
                 is_xahuman_split  => $contains_human ? 0 : $contains_xahuman,
                 is_yhuman_split  => $contains_human ? 0 : $contains_yhuman,
                 human_split_type =>  $contains_human ? 'human' : ( $contains_xahuman ? 'xahuman' : ( $contains_yhuman ? 'yhuman' : undef ) ),
                 is_spiked_phix  => $lims->spiked_phix_tag_index ? 1 : 0,
                 spiked_phix_tag_index  => $spiked_phix_tag_index,
               };

  foreach my $meta_name ( @LIMS_PROPERTIES ){
      if ($meta_name eq q{sample_name} || $meta_name eq q{study_name} ){
        next;
      }
      $names->{$meta_name} = [sort keys %{$data->{$meta_name}}];
  }

  return $names;
}

=head2 bam_fully_archived_for_deletion

simply checking bam file in archive for deletion

=cut

sub bam_fully_archived_for_deletion { ## no critic (Subroutines::ProhibitExcessComplexit)
    my $self = shift;

    my $fully_archived = 1;

    # get a list of archived files
    my $bam_list_archived;
    eval{
        $bam_list_archived = npg_common::irods::Loader->new(file=>'none')
                            ->get_collection_file_list($self->collection());
	     1;
    } or do {
        carp "Problems to check irods: $EVAL_ERROR";
	     return 0;
    };

    my @file_list = keys %{$bam_list_archived};
    my %bam_list = map {$_=> 1 } grep { /bam$/mxs } @file_list;
    my %bai_list = map {$_=> 1 } grep { /bai$/mxs } @file_list;

    # check bai files match bam files
    if(! $self->_check_bam_index_file(\%bam_list, \%bai_list)){
       $fully_archived = 0;
    }

    #check no bam file missing in irods
    my %bam_list_by_lane = ();
    my %bam_list_by_lane_split = ();

    my %bam_list_by_lane_tag = ();
    my %bam_list_by_lane_split_tag = ();

    my %bam_list_by_lane_spiked_phix = ();
    my %bam_list_by_lane_spiked_phix_tag = ();

    foreach my $file_name (keys %bam_list){

        my ($lane)        = $file_name =~ m/ \d+ _ (\d) /mxs;
        my ($tag_index)   = $file_name =~ m/ [#] (\d+) /mxs;
        my ($split) = $file_name =~ m/ _ (yhuman|xahuman|human|phix) /mxs;

        if(!$lane){
           croak "There is no lane number in bam file $file_name";
        }

        if ( !$split ){
           if (!defined $tag_index){
             $bam_list_by_lane{$lane}++;
           } else{
             $bam_list_by_lane_tag{$lane}->{$tag_index}++;
           }
        } elsif ( $split =~ /human/mxs ) {
           if (!defined $tag_index){
             $bam_list_by_lane_split{$lane}->{human_split}++;
           } else{
             $bam_list_by_lane_split_tag{$lane}->{human_split}->{$tag_index}++;
           }
        }elsif ( $split eq q{phix} ){
           if (!defined $tag_index){
             $bam_list_by_lane_spiked_phix{$lane}++;
           } else{
             $bam_list_by_lane_spiked_phix_tag{$lane}->{$tag_index}++;
           }
        }
    }

    my $positions = $self->positions();

    if( ! $self->positions() || ( scalar @{$self->positions()} == 0) ){
       $positions = [1..$DEFAULT_LANE_NUMBERS];
    }

    foreach my $position (@{$positions}){

          my $st_info_lane = $self->get_study_library_sample_names($self->id_run(), $position);
          my $is_phix_control_lane = $st_info_lane ->{is_phix_control};
          my $is_human_split_lane = $st_info_lane->{is_human_split} ? 1 : $st_info_lane->{is_xahuman_split} ? 1 : $st_info_lane->{is_yhuman_split} ? 1 : 0;
          my $is_spiked_phix = $st_info_lane->{is_spiked_phix};
          my $spiked_phix_tag_index = $st_info_lane->{spiked_phix_tag_index};
          my $is_multiplexed = $self->is_multiplexed_lane($position);

          if($is_phix_control_lane){
             next;
          }

          if( $is_spiked_phix && !$self->alignment_filter_in_plex()  && ! $bam_list_by_lane_spiked_phix{$position} ){
              $self->log("Lane $position spiked phix bam file missing");
              $fully_archived = 0;
          }

          if(!$is_multiplexed){
             my $number_bam_files = $bam_list_by_lane{$position};
             my $number_bam_files_split = $bam_list_by_lane_split{$position}->{human_split};

             if( ( $is_human_split_lane && $number_bam_files_split !=1 )
                || (!$is_human_split_lane &&  !$number_bam_files ) ){

                   $self->log("Lane $position bam file missing");
                   $fully_archived = 0;
              }
          }else {
              my $tag_index_list = $self->get_tag_index_list($position);

              foreach my $tag_index (@{$tag_index_list}){

                 if( $spiked_phix_tag_index
                     && $self->alignment_filter_in_plex()
                     && $tag_index != $spiked_phix_tag_index
                     && ! $bam_list_by_lane_spiked_phix_tag{$position}->{$tag_index} ){
                   $self->log("Lane $position tag $tag_index spiked phix bam file missing");
                   $fully_archived = 0;
                 }

                 if( !$bam_list_by_lane_tag{$position}->{$tag_index}
                     && !$bam_list_by_lane_split_tag{$position}->{human_split}->{$tag_index} ){
                         $self->log("Lane $position tag $tag_index bam file missing");
                         $fully_archived = 0;
                 }
              }

          }
    }

    $self->log('The bam list in irods is correct against lims system.');

    # check archived bam files with correct md5
    my $bam_list_to_archive = $self->file_list();
    if( ! $self->check_bam_list_md5(\%bam_list, $bam_list_to_archive) ){
        $self->log('The number or md5 values of bam files on irods is not the same as on staging');
        return 0;
    }

    return $fully_archived;
}

=head2 check_bam_list_md5

compare md5 values between two bam list from staging and irods

=cut

sub check_bam_list_md5 {
   my ( $self, $irods_bam_list, $staging_bam_list ) = @_;

   my $number_bam_on_irods = keys %{$irods_bam_list};
   my $number_bam_on_staging = scalar @{ $staging_bam_list};
   if( $number_bam_on_irods != $number_bam_on_staging ){
      $self->log("Number of bam files is different. iRODs: $number_bam_on_irods, Staging: $number_bam_on_staging");
      return;
   }

   my $md5_list_irods   = $self->_get_irods_bam_list_md5s( $irods_bam_list );
   my $md5_list_staging = $self->_get_staging_bam_list_md5s( $staging_bam_list );
   my $md5_correct = 1;
   foreach my $bam (keys %{$md5_list_irods}){
      my $md5_irods = $md5_list_irods->{$bam};
      my $md5_staging = $md5_list_staging->{$bam};
      if( !$md5_irods || !$md5_staging || $md5_irods ne $md5_staging){
         $self->log("BAM $bam md5 wrong: $md5_irods not match $md5_staging");
         $md5_correct = 0;
      }
   }
   return $md5_correct;
}

sub _get_irods_bam_list_md5s {
   my ($self, $irods_bam_list) = @_;

   my $md5_list = {};

   foreach my $irods_bam ( keys %{$irods_bam_list} ) {
      my $file_name =  basename($irods_bam);
      my $loader = npg_common::irods::Loader->new();
      my $md5 = $loader->_get_irods_md5( $self->collection() . q{/} . $irods_bam );
      $md5_list->{$file_name} = $md5;
   }
   return $md5_list;
}

sub _get_staging_bam_list_md5s {
   my ($self, $staging_bam_list) = @_;

   my $md5_list = {};

   foreach my $bam ( @{$staging_bam_list} ) {
      my $file_name =  basename($bam);
      my $loader = npg_common::irods::Loader->new();
      my $md5 = $loader->_get_file_md5($bam);
      $md5_list->{$file_name} = $md5;
   }
   return $md5_list;
}


sub _check_bam_index_file {
    my ($self, $bam_list_ref, $bai_list_ref) = @_;

    if($self->exclude_bam_index()){
       return 1;
    }

    my $index_file_found = 1;
    foreach my $bam ( keys %{$bam_list_ref} ){

       my $total_reads = $self->get_number_of_reads_irods_meta($bam);
       if(( defined $total_reads && ($total_reads == 0) ) || $self->_bam_no_alignment_in_irods($bam) ){
          next;
       }

       my $bai = $bam;
       $bai =~ s/bam$/bai/mxs;
       if(!$bai_list_ref->{$bai}){
          $self->log("Index file for $bam not exist");
          $index_file_found = 0;
       }
    }

    return $index_file_found;
}

=head2 bam_fully_archived

check every bam file fully archived

=cut

sub bam_fully_archived {
    my $self = shift;

    my $fully_archived = 1;

    my $bam_list_archived = npg_common::irods::Loader->new(file=>'none')
                            ->get_collection_file_list($self->collection());
    my $bam_list_to_archive = $self-> _get_bam_list_to_archive();

    foreach my $bam (@{$bam_list_to_archive}){
       if(!$bam_list_archived->{$bam}){
          $fully_archived = 0;
          $self->log("File Missing in iRods: $bam");
       }

       if( $self->exclude_bam_index() ){
           next;
       }

       if( ! $self->get_number_of_reads($bam) ){
           next;
       }

       my $bam_index_file = $bam;
       $bam_index_file =~ s/bam$/bai/mxs;
       if(!$bam_list_archived->{$bam_index_file} && !$self->_bam_no_alignment_in_irods( $bam )){
          $fully_archived = 0;
          $self->log("File Missing in iRods: $bam_index_file");
       }

    }

    return $fully_archived;
}

sub _get_bam_list_to_archive {
    my $self = shift;

    $self->log('Getting a list of bam files which should be archived');

    my @bam_list = ();

    my $positions = $self->positions();
    if ( !$self->positions() || ( scalar @{$self->positions()} == 0) ) {
       $positions = [1..$DEFAULT_LANE_NUMBERS];
    }

    foreach my $position (@{$positions}){
       push @bam_list, @{$self->_bam_file_list_per_lane($position)};
    }
    return \@bam_list;
}

sub _bam_file_list_per_lane {
   my ($self, $position) = @_;

   my @bam_list_per_lane = ();

   my $st_info_lane = $self->get_study_library_sample_names($self->id_run(), $position);
   my $is_phix_control_lane = $st_info_lane ->{is_phix_control};

   my $is_spiked_phix      = $st_info_lane->{is_spiked_phix};
   my $spiked_phix_tag_index = $st_info_lane->{spiked_phix_tag_index};

   my $split = $st_info_lane->{'human_split_type'};

   if ( !$self->is_multiplexed_lane($position) && !$is_phix_control_lane ) {
          if( $split){
             push @bam_list_per_lane, $self->_bam_file_name({id_run      => $self->id_run(),
                                                             position    => $position,
                                                             split       => $split});
          }
          push @bam_list_per_lane, $self->_bam_file_name({id_run   => $self->id_run(),
                                                             position => $position});
    } elsif (!$is_phix_control_lane){

          my $tag_index_list = $self->get_tag_index_list($position);

          foreach my $tag_index (@{$tag_index_list}){
             push @bam_list_per_lane, $self->_bam_file_name({id_run   => $self->id_run(),
                                                                position => $position,
                                                                tag_index=> $tag_index,});
             if ( !$spiked_phix_tag_index || $tag_index != $spiked_phix_tag_index) {
                if ( $split ) {
                   push @bam_list_per_lane, $self->_bam_file_name({id_run      => $self->id_run(),
                                                                position    => $position,
                                                                split       => $split,
                                                                tag_index   => $tag_index,});
                }
                if ( $spiked_phix_tag_index && $self->alignment_filter_in_plex() ) {
                   push @bam_list_per_lane, $self->_bam_file_name({
                                                        id_run => $self->id_run(),
                                                        position    => $position,
                                                        tag_index=> $tag_index,
                                                        split => 'phix',});
                }
	     }

          }
    }

    if( $is_spiked_phix && !$self->alignment_filter_in_plex() ) {
        push @bam_list_per_lane, $self->_bam_file_name({id_run      => $self->id_run(),
                                                        position    => $position,
                                                        split => 'phix'});
    }

    return \@bam_list_per_lane;
}

sub _bam_file_name {
   my ($self, $args_ref) = @_;

   my $id_run      = $args_ref->{id_run};
   my $position    = $args_ref->{position};
   my $split       = $args_ref->{split};
   my $tag_index   = $args_ref->{tag_index};

   my $bam_file_name = $id_run.q{_}.$position;

   if( $self->alignment_filter_in_plex() ){
       if( defined $tag_index ){
           $bam_file_name .= q{#}.$tag_index;
       }
       if( $split ){
           $bam_file_name .= q{_}.$split;
       }
   }else{
      if( $split ){
         $bam_file_name .= q{_}.$split;
      }
      if( defined $tag_index ){
         $bam_file_name .= q{#}.$tag_index;
      }
   }

   return $bam_file_name.q{.bam};
}

sub _get_ref_from_bwa_pg {
   my ($self, $bwa_pg_line) = @_;

   my @ref = grep {/\/(nfs|lustre)\/\S+\/references/mxs} split /\s/mxs, $bwa_pg_line;
   if(scalar @ref){
      return $ref[0];
   }
   return;
}

=head2 get_bam_reference

read bam header using samtools and check bwa or bwa_aln PG to get the reference used

=cut

sub get_bam_reference {
   my ( $self, $bam ) = @_;

   my $samtools_view_cmd = $self->samtools_cmd . qq[ view -H $bam];
   $self->log($samtools_view_cmd);

   my $pid = open3( undef, my $sam_header_fh, undef, $samtools_view_cmd);

   my $last_bwa_pg_line;
   my $sq_tag_found = 0;
   while( my $line = <$sam_header_fh> ) {
      if( $line =~ /^\@PG/mxs && $line =~/ID\:bwa/mxs && $line !~ /ID\:bwa_sam/mxs && $line =~ /\tCL\:/mxs) {
          chomp $line;
          $last_bwa_pg_line = $line;
      }elsif($line =~ /^\@SQ\t/mxs){
          $sq_tag_found = 1;
      }
   }

   my $reference;
   if($last_bwa_pg_line){
      $reference = $self->_get_ref_from_bwa_pg( $last_bwa_pg_line );
   }
   waitpid $pid, 0;

   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak qq{Failed: $samtools_view_cmd};
   }

   close $sam_header_fh or croak "can not close samtools view output: $ERRNO";

   if($reference){
     $self->log("Reference found from bam header: $reference");
   }else{
       $self->log('NO reference found from bam header.');
   }
   return {reference => $reference, alignment => $sq_tag_found };
}

=head2 get_number_of_reads

get number of reads of the given bam file from cache, bam flagstats json or running samtools

=cut

sub get_number_of_reads {
   my ($self, $bam) = @_;

   my ($bam_name, $path, $suffix) = fileparse($bam, q{.bam});

   if( defined $self->number_of_reads_list()->{$bam_name} ){
      return $self->number_of_reads_list()->{$bam_name};
   }

   my $total_reads = $self->get_number_of_reads_qc( $bam );

   if( ! defined $total_reads ){
      $total_reads = $self->get_number_of_reads_samtools( $bam );
   }

   if(defined $total_reads) {

      $self->number_of_reads_list()->{$bam_name} = $total_reads;
   }

   return $total_reads;
}

=head2 get_number_of_reads_qc

get number of reads of the given bam file bam flagstats json

=cut

sub get_number_of_reads_qc {
   my ($self, $bam) = @_;

   my $num_total_reads;

   my ($name, $path, $suffix) = fileparse($bam, q{.bam});

   my $bam_flagstats_file = File::Spec->catfile ( $path, q{qc}, $name.q{_bam_flagstats.json} );

   if(-e $bam_flagstats_file ){

      my $flag_stats = from_json(slurp $bam_flagstats_file);
      $num_total_reads = $flag_stats->{num_total_reads};
   }

   return $num_total_reads;
}

=head2 get_number_of_reads_samtools

get number of reads of the given bam file from running samtools

=cut

sub get_number_of_reads_samtools {
   my ( $self, $bam ) = @_;

   my $samtools_flagstat_cmd = $self->samtools_cmd . qq[ flagstat $bam];
   $self->log($samtools_flagstat_cmd);

   my $pid = open3( undef, my $sam_flagstat_fh, undef, $samtools_flagstat_cmd);

   my $total_reads = 0;

   while( my $line = <$sam_flagstat_fh> ) {

      chomp $line;
      if( $line =~ /^(\d+)[ ].*in[ ]total/mxs ){
         $total_reads = $1;
         last;
      }
   }

   waitpid $pid, 0;

   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak qq{Failed: $samtools_flagstat_cmd};
   }

   close $sam_flagstat_fh or croak "can not close samtools flagstat output: $ERRNO";

   return $total_reads;
}

=head2 get_number_of_reads_irods_meta

get number of reads of the given bam file from irods meta data

=cut

sub get_number_of_reads_irods_meta {
   my ( $self, $bam ) = @_;

   if( exists $self->number_of_reads_list_irods()->{$bam} ){
       return $self->number_of_reads_list_irods()->{$bam};
   }

   my ($id_run) = $bam =~ /^(\d+)_/mxs;
   $id_run += 0;

   if(!$id_run){
      croak "no id_run from bam file name: $bam";
   }

   my $irods_bam = File::Spec->catfile ( $self->collection(), $bam);

   my $loader = npg_common::irods::Loader->new(file=>'none');
   my $meta_list = $loader->_check_meta_data ($irods_bam);

   $self->_bam_indexed_cache_irods()->{$bam} = (exists($meta_list->{alignment}) && exists($meta_list->{alignment}->{1})) ? 1 : 0;

   if(  !$meta_list->{total_reads} ){
      $self->number_of_reads_list_irods()->{$bam} = undef;
      return;
   }

   my @total_reads_values = keys %{ $meta_list->{total_reads} };

   if(scalar @total_reads_values > 1) {
     croak "multiple total reads value in irods for bam: $bam";
   }elsif(scalar @total_reads_values){
     $self->number_of_reads_list_irods()->{$bam} = $total_reads_values[0];
     return $total_reads_values[0];
   }

   return;
}

sub _bam_no_alignment {
   my ($self, $bam_file) = @_;

   my $bam_indexed_cache =  $self->_bam_indexed_cache();

   if( exists $bam_indexed_cache->{$bam_file} ){
      return $bam_indexed_cache->{$bam_file};
   }

   my $mark_duplicates_obj = $self->_bam_to_mark_duplicate($bam_file);

   $bam_indexed_cache->{$bam_file} = $mark_duplicates_obj->no_alignment();

   return $bam_indexed_cache->{$bam_file};
}

sub _bam_to_mark_duplicate {
   my ($self, $bam_file) = @_;

   return npg_common::sequence::BAM_MarkDuplicate->new(
                  input_bam     => $bam_file,
                  output_bam    => q{no_output},
                  metrics_json  => q{no_json_output},
                  $self->resolved_paths,
              );
}

sub _bam_no_alignment_in_irods {
   my ($self, $bam) = @_;

   my $bam_indexed_cache_irods =  $self->_bam_indexed_cache_irods();

   if( exists $bam_indexed_cache_irods->{$bam} ){
      return !$bam_indexed_cache_irods->{$bam};
   }

   $self->get_number_of_reads_irods_meta($bam);

   return !$bam_indexed_cache_irods->{$bam};
}

no Moose;

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item MooseX::Getopt

=item Carp

=item English -no_match_vars

=item Readonly

=item File::Basename

=item File::Spec

=item st::api::lims

=item npg_common::irods::Loader

=item npg_tracking::illumina::run::short_info

=item npg_tracking::illumina::run::folder

=item npg_common::roles::log

=item npg_tracking::illumina::run::long_info

=item npg_common::roles::run::lane::tag_info

=item npg_common::roles::software_location

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by Guoying Qi

This file is part of NPG.

NPG is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
