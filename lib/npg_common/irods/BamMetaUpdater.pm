#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-07-01
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::BamMetaUpdater;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use IPC::Open3;
use POSIX qw(strftime);
use List::MoreUtils qw{any};

use npg_common::irods::Loader;
use npg_warehouse::Schema;

with qw /
         MooseX::Getopt
         npg_common::roles::software_location
         npg_common::roles::log
        /;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $EXIT_CODE_SHIFT    => 8;

Readonly::Scalar our $NO_METADATA_FILE   => q{no_lims_data};
Readonly::Scalar my  $NOT_SPIKED         => -1;
Readonly our %IGNORE_RUN_LIST => ();
Readonly our %IGNORE_STUDY_LIST  => ();
Readonly our %IGNORE_SAMPLE_LIST => ();


## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::BamMetaUpdater

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 SUBROUTINES/METHODS

=head2 id_run

id_run

=cut
has 'id_run'       => ( isa           => 'Int',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'the id_run filter',
                      );

=head2 lane

lane

=cut
has 'lane'         => ( isa           => 'Int',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'the lane number filter',
                      );


=head2 tag_index

tag_index

=cut
has 'tag_index'    => ( isa           => 'Int',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'the tag index number filter',
                      );


=head2 min_id_run

min_id_run

=cut

has 'min_id_run'   => ( isa           => 'Int',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'the minimum id_run filter',
                      );

=head2 max_id_run

max_id_run

=cut

has 'max_id_run'   => ( isa           => 'Int',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'the maximum id_run filter',
                      );

=head2 check_bam

check bam header for meta data change

=cut
has 'check_bam'=>        ( isa           => 'Bool',
                           is            => 'ro',
                           required      => 0,
                           documentation => 'check bam header for meta data change',
                         );

=head2 check_md5

check md5 value in irods meta data match the value reported by ils -L command 

=cut
has 'check_md5'=>  ( isa           => 'Bool',
                     is            => 'ro',
                     required      => 0,
                     documentation => 'check md5 value in irods meta data match the value reported by ils command ',
                   );

=head2 not_check_warehouse

do not check warehouse for lims data change

=cut
has 'not_check_warehouse'=>  ( isa           => 'Bool',
                               is            => 'ro',
                               required      => 0,
                               documentation => 'do not check warehouse for lims data change',
                   );

=head2 verbose

Verbose flag

=cut
has 'verbose'      => ( isa        => 'Bool',
                        is         => 'ro',
                        required   => 0,
                        default    => 0,
                      );

=head2 not_dry_run

not_dry_run

=cut

has 'not_dry_run'  => ( isa           => 'Bool',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'imeta sub commands will be run and irods meta data will be updated',
                      );

=head2 output

output

=cut
has 'output'        => ( isa          => 'Str',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'The file name to store imeta sub commands which will be printed out if not given',
                      );

=head2 _schema_wh

DBIx schema object for the warehouse database

=cut
has '_schema_wh'  =>  ( isa        => 'npg_warehouse::Schema',
                        is         => 'ro',
                        required   => 0,
                        lazy_build => 1,
                      );

sub _build__schema_wh {
    my $self = shift;
    my $schema = npg_warehouse::Schema->connect();
    if($self->verbose) {
        $self->log( q[Connected to the warehouse db, schema object ] . $schema );
    }
    return $schema;
}

=head2 process

main method to call, get each bam file, compare meta data in irods with warehouse

update them if different

=cut

sub process{ ## no critic ( Subroutines::ProhibitExcessComplexity )
  my $self = shift;

  my $output_fh;
  if($self->output()){
     open $output_fh, q{>}, $self->output() or croak q{can not open the output file }. $self->output();## no critic (InputOutput::RequireBriefOpen )
  }

   my  $file_list = $self->get_file_list();

   $self->verbose() && $self->log('The number files in irods to check: ' . @{$file_list});

   my $irods_loader = npg_common::irods::Loader->new();
   my @imeta_commands_all;

   foreach my $irods_file (@{$file_list}){

           $self->verbose() && $self->log( $irods_file);

           my ($id_run, $position, $tag_index, $alignment_filter)
                           = $self->parsing_bam_filename($irods_file);
	   if (!$id_run) {
	      carp "Failed to infer id_run from the file name '$irods_file'; skipping to next file";
	      next;
	   }

           #ignore some runs
           if( $IGNORE_RUN_LIST{$id_run} ){
              $self->verbose && $self->log("Run $id_run ignored ");
              next;
           }

           my $irods_meta;
           eval{
              $irods_meta = $self->get_current_irods_meta($irods_file);
              1;
           } or do {
              carp $EVAL_ERROR;
              next;
           };

           if ($irods_meta->{$NO_METADATA_FILE}) { # this flag is set manually
              $self->verbose && $self->log("File $irods_file ignored since it has '$NO_METADATA_FILE' flag set");
	      next;
	   }

           my $meta_list_warehouse = {};
           if(!$self->not_check_warehouse()){
             $meta_list_warehouse = $self->get_warehouse_npg_information($id_run, $position, $tag_index);
           }

           if( $self->check_bam() ){
              my $meta_list_from_bam  = $self->check_bam_header($irods_file);
              $meta_list_warehouse = { %{$meta_list_warehouse}, %{$meta_list_from_bam} };
           }

           if( $self->check_md5() ){
              my $md5_from_irods   = $self->get_irods_file_md5($irods_file);
              $meta_list_warehouse = { %{$meta_list_warehouse}, %{$md5_from_irods } };
           }

           my $total_reads;
           my $total_reads_values_hashref = $irods_meta->{total_reads};
           if($total_reads_values_hashref){
              ($total_reads) = keys %{$total_reads_values_hashref};
           }
           if(defined $total_reads && $total_reads == 0 && exists $meta_list_warehouse->{alignment} ){
               $meta_list_warehouse->{alignment} = [0];
           }

           #ignore some cases where some rows should be deleted from warehouse npg_plex_information table
           if (defined $tag_index
                && $tag_index == 0
                && $meta_list_warehouse->{sample}
                && ( scalar keys %{$irods_meta->{sample}} < scalar @{$meta_list_warehouse->{sample} } )
              ){
              $self->verbose && $self->log('Some old records in warehouse should be deleted');
              next;
           }

           if(!$meta_list_warehouse){
              $self->verbose && $self->log("THERE IS NO INFORMATION in warehouse for $irods_file");
              next;
           }

           my $imeta_commands = $self->check_meta_data($irods_file, $irods_meta, $meta_list_warehouse);

           if( scalar @{$imeta_commands} ){

              $self->log("NOT OK: $irods_file");

              if($self->output()){
                 print {$output_fh} (join "\n", @{$imeta_commands})."\n" or croak 'cannot print to output';
              }else{
                 push @imeta_commands_all, @{$imeta_commands};
              }

              if($self->not_dry_run()){
                my $imeta_commands_string = join qq{\n}, @{$imeta_commands};
                $self->log("imeta commands to run:\n$imeta_commands_string");
                $irods_loader->_run_imeta_commands( $imeta_commands_string );
                $self->log( qq{irods meta data updated for $irods_file!\n} );
              }
           }

   }

   if(!$self->not_dry_run()){
      $self->log('irods meta data not updated, see output file to check anything changed');
   }

   if($output_fh){
      close $output_fh or croak 'can not close output';
   }elsif(!$self->not_dry_run()){
      my $imeta_commands_string_all = join qq{\n}, @imeta_commands_all;
      $self->log("\n". $imeta_commands_string_all);
   }

   return 1;
}

=head2 get_file_list

method to query irods to get a list of bam files

=cut
sub get_file_list{
  my $self = shift;

  my $imeta_qu_bam_list_cmd = qq{imeta qu -z $npg_common::irods::Loader::DEFAULT_ZONE -d type = bam};

  if($self->id_run()){
     $imeta_qu_bam_list_cmd .= q{ and id_run = } . $self->id_run();
  }
  if($self->lane()){
     $imeta_qu_bam_list_cmd .= q{ and lane = } . $self->lane();
  }
  if(defined $self->tag_index()){
     $imeta_qu_bam_list_cmd .= q{ and tag_index = } . $self->tag_index();
  }
  if( $self->min_id_run() ) {
     $imeta_qu_bam_list_cmd .= q{ and id_run 'n>=' } . $self->min_id_run();
  }
  if( $self->max_id_run() ) {
     $imeta_qu_bam_list_cmd .= q{ and id_run 'n<=' } . +$self->max_id_run();
  }

  $self->verbose() && $self->log($imeta_qu_bam_list_cmd);

  my @file_list = ();

  my $pid = open3( undef, my $imeta_qu_out_fh, undef, $imeta_qu_bam_list_cmd);

  my ($collection, $file);

  while (my $line = <$imeta_qu_out_fh> ) {

     chomp $line;

     if($line =~ /^[-]{4}/mxs){

        $collection = undef;
        $file = undef;
     }elsif($line =~ /collection[:][ ](.+)/mxs){

        $collection = $1;
     }elsif ($line =~ /dataObj[:][ ](.+)/mxs){

        $file = $1;
        if(! $collection || !$file){
           croak "Problems with irods output from $imeta_qu_bam_list_cmd";
        } else{

           my $irods_file = $collection.q{/}.$file;

           push @file_list, $irods_file;

        }
     }
   }

   waitpid $pid, 0;

   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak "Failed: $imeta_qu_bam_list_cmd";
   }

   close $imeta_qu_out_fh or croak "can not close imeta command output: $ERRNO";

   return \@file_list;
}

=head2 check_meta_data

compare current irods meta data for one file with warehouse value, return a list of imeta commands if different

=cut
sub check_meta_data {
   my ($self, $irods_file, $irods_meta, $meta_list_warehouse) = @_;

   my @imeta_commands = ();

   foreach my $meta ( keys %{$meta_list_warehouse} ){

       my $error_messages;

       my @values_warehouse_not_trimmed = @{$meta_list_warehouse->{$meta}};
       my @values_warehouse = ();
       foreach my $value (@values_warehouse_not_trimmed){
          chomp $value;
          if(defined $value && $value ne q{}){
             $value =~ s/^\s+//mxs;
             $value =~ s/\s+$//mxs;
             push @values_warehouse, $value;
          }
       }

       my $values_in_irods = $irods_meta->{$meta};

       if(scalar @values_warehouse == 0){
          $self->verbose && $self->log("THERE IS NO $meta INFORMATION in warehouse");
          next;
       }

       if(scalar @values_warehouse != scalar keys %{$values_in_irods} ){
           $error_messages = "$meta: The number of values different. WareHouse: "
                                 . (join q{,}, @values_warehouse )
                                 . q{. irods: }
                                 . (join q{,}, keys %{$values_in_irods});
           $self->log( $error_messages );
           push @imeta_commands, @{$self->change_meta_data($irods_file, $meta,$values_in_irods, \@values_warehouse)};
           next;
       }

       my $matched = 1;
       foreach my $value (@values_warehouse){

          chomp $value;
          #we change the value when we adding meta data to irods at the begining
          if($value =~ /"/mxs && $value =~ /'/mxs ){
             $value =~ s/"/'/gmxs;
          }

          if(! $values_in_irods->{$value}){
             $error_messages = qq{$meta: $value missing in irods. irods: } .(join q{,}, keys %{$values_in_irods});
             $self->verbose && $self->log( $error_messages );
             $matched = 0;
          }
       }

       if(!$matched){

          my @meta_history_values = keys %{$irods_meta->{$meta . q{_history}}};

          my @current_meta_values = keys %{$values_in_irods};
          @current_meta_values = sort @current_meta_values;
          my $current_meta_values_string = join q{,}, @current_meta_values;

          my %existd_history = map{ $self->_remove_history_date($_)=>1 } @meta_history_values;

          my $no_add_history = 0;
          if ($existd_history{$current_meta_values_string} ){

             $self->log("$meta meta values are different but history record already added: " . join "\n", keys %existd_history );
             $no_add_history = 1;
          }

          push @imeta_commands, @{$self->change_meta_data($irods_file, $meta, $values_in_irods, \@values_warehouse, $no_add_history)};
       }
   }

   return \@imeta_commands;
}

sub _remove_history_date {
    my ($self, $history_value) = @_;

    $history_value =~ s/\[\d{4}\-\d{2}\-\d{2}T\d{2}\:\d{2}\:\d{2}\]\s+//mxs;

    return $history_value;
}

=head2 change_meta_data

to form the imeta commands

=cut
sub change_meta_data { ## no critic (Subroutines::ProhibitManyArgs)
   my ($self, $irods_file, $meta, $old_values_hash, $new_values_array, $no_add_history) = @_;

   my @imeta_commands = ();

   my @old_values = keys %{$old_values_hash};

   @old_values = sort @old_values;

   if(scalar @old_values ){

      push @imeta_commands, qq{rmw -d $irods_file $meta "%"};

      if(!$no_add_history){

         my $old_values_string = (strftime q{[%Y-%m-%dT%H:%M:%S] }, localtime). (join q{,}, @old_values);
         my $meta_name = $meta . q{_history};
         push @imeta_commands, _imeta_add_command($irods_file, $meta_name, $old_values_string);
      }
   }

   if(scalar @{$new_values_array}){
      my @new_add_command = map { _imeta_add_command($irods_file, $meta, $_)  } @{$new_values_array};
      push @imeta_commands, @new_add_command;
   }

   return \@imeta_commands;
}

sub _imeta_add_command {
   my ($irods_file, $meta_name, $meta_value) = @_;

   chomp $meta_value;

   my $quoted = $meta_value;
   if($quoted =~ /"/mxs && $quoted !~ /'/mxs){
    $quoted = qq{'$quoted'};
  }elsif($quoted =~ /'/mxs && $quoted !~ /"/mxs){
    $quoted = qq{"$quoted"};
  }else{
    $quoted =~ s/"/'/gmxs;
    $quoted = qq{"$quoted"};
  }

  return qq{add -d $irods_file $meta_name $quoted};
}

=head2 get_warehouse_npg_information 

get warehouse information

=cut

sub get_warehouse_npg_information {## no critic ( Subroutines::ProhibitExcessComplexity )
   my ($self, $id_run, $position, $tag_index) = @_;

   my $schema = $self->_schema_wh();
   my $npg_information = $schema->resultset('NpgInformation')->find(
        {
           id_run    => $id_run,
           position  => $position,
        },
        { key => 'id_run_position' }
   );

   if (!$npg_information) {
     croak "Failed to retrieve a row from npg_information table: id_run $id_run, position $position";
   }

   my $is_paired_read = $npg_information->paired_read();
   my $spiked_phix_phex_index = $npg_information->spike_tag_index();
   $spiked_phix_phex_index ||= $NOT_SPIKED;

   my $plex_query = {
           id_run    => $id_run,
           position  => $position,
                    };

   if ($tag_index) {
      $plex_query->{tag_index} = $tag_index;
   }

   my @npg_plex_information = $schema->resultset('NpgPlexInformation')->search( $plex_query );


   if( !$npg_information && scalar @npg_plex_information == 0 ) {
      return;
   }


   my %samples = ();
   my %libraries = ();
   my %studies =();
   my %library_ids = ();
   my %sample_public_names = ();
   my %sample_accession_numbers = ();
   my %sample_ids = ();
   my %sample_common_names = ();
   my %study_ids = ();
   my %study_titles = ();
   my %study_accession_numbers = ();

   my @manual_qc = ();
   my $manual_qc_state = $npg_information->manual_qc();
   if( defined $manual_qc_state ) {
      push @manual_qc, $manual_qc_state;
   }

   #study information
   my @study_ids;
   if( $tag_index && $tag_index == $spiked_phix_phex_index ){
       @study_ids = map { $_->study_id() } @npg_plex_information;
   }else{
       @study_ids = map { $_->study_id() } grep {$_->tag_index != $spiked_phix_phex_index} @npg_plex_information;
   }

   if(scalar @study_ids == 0 ){
       push @study_ids, $npg_information->study_id();
   }

   foreach my $study_id ( @study_ids ){

       if( !$study_id ||  $IGNORE_STUDY_LIST{ $study_id } ){
         next;
       }

       $study_ids{$study_id}++;

       my $studies_rs = $schema->resultset('CurrentStudy')->search( {internal_id=>$study_id} );
       my $num_studies = $studies_rs->count;
       if (!$num_studies) {
	 croak "No record for study $study_id in the CurrentStudy table of the warehouse";
       }
       if ($num_studies > 1) {
	 croak "Multiple records for study $study_id in the CurrentStudy table of the warehouse";
       }
       my $study = $studies_rs->next;
       my $name = $study->name();
       if($name){
          $studies{$name}++;
       }

       my $study_title = $study->study_title();
       if($study_title){
          $study_titles{$study_title}++;
       }

       my $accession_number = $study->accession_number();
       if($accession_number){
          $study_accession_numbers{$accession_number}++;
       }

   }

   #library 
   if(! $tag_index ){
      my $library = $npg_information->asset_name();
      my $library_id = $npg_information->asset_id();
      $library ||= $library_id;
      if($library){
         $libraries{$library}++;
      }
      if($library_id){
         $library_ids{$library_id}++;
      }
   }else{
      foreach my $npg_info (@npg_plex_information){
        my $library = $npg_info->asset_name();
        my $library_id = $npg_info->asset_id();
        $library ||= $library_id;
        if($library){
           $libraries{$library}++;
        }
        if($library_id){
           $library_ids{$library_id}++;
        }
      }
   }

   #sample names
   my @sample_objects = ();
   if( ! defined $tag_index ){
        my $sample = $npg_information->sample();
        if($sample){
           push @sample_objects, $sample;
        }
   }

   ## sample names for plex or pool lane
   if(defined $tag_index || scalar @sample_objects == 0){

     foreach my $npg_info (@npg_plex_information){

        my $sample = $npg_info->sample();
        my $plex_tag_index = $npg_info->tag_index();

        if( $plex_tag_index == $spiked_phix_phex_index
           && ( ( $tag_index && $tag_index != $spiked_phix_phex_index ) || !$tag_index ) ) {
           next;
        }

        if($sample){
           push @sample_objects, $sample;
        }

     }
   }

   my $consent_withdrawn = 0;
   foreach my $sample (@sample_objects){

           my $sample_id = $sample->internal_id();
           if( $IGNORE_SAMPLE_LIST{$sample_id} ){
             next;
           }

           $sample_ids{$sample_id}++;

           my $sample_name = $sample->name();

           if($sample_name) {
              $samples{$sample_name}++;
           }

           my $public_name = $sample->public_name();

           if($public_name){
              $sample_public_names{$public_name}++;
           }

           my $accession_number = $sample->accession_number();

           if($accession_number ){
              $sample_accession_numbers{$accession_number}++;
           }

           my $common_name = $sample->common_name();

           if($common_name){
	      chomp $common_name;
	      $common_name =~ s/^\s+//mxs;
	      $common_name =~ s/\s+$//mxs;
              $sample_common_names{$common_name}++;
           }
           if( !$consent_withdrawn && $sample->consent_withdrawn ) {
	      $consent_withdrawn = 1;
	   }
   }

   my $meta_data_list = {sample  => [keys %samples],
                         library => [keys %libraries],
                         study   => [keys %studies],
                         library_id => [keys %library_ids],
                         sample_public_name => [keys %sample_public_names],
                         sample_accession_number => [keys %sample_accession_numbers],
                         sample_public_name => [keys %sample_public_names],
                         sample_common_name => [keys %sample_common_names],
                         sample_id          => [keys %sample_ids],
                         study_id           => [keys %study_ids],
                         study_title        => [keys %study_titles],
                         study_accession_number => [keys %study_accession_numbers],
                         manual_qc          => \@manual_qc,
                         is_paired_read     => [$is_paired_read],
                        };
   if ($consent_withdrawn) {
     $meta_data_list->{'sample_consent_withdrawn'} = [1];
     $meta_data_list->{'target'} = [0];
   }

   return $meta_data_list;
}

=head2 get_current_irods_meta

get current irods meta data

=cut
sub get_current_irods_meta {
   my ($self, $irods_file) = @_;
   return  npg_common::irods::Loader->new()->_check_meta_data($irods_file);
}

=head2 parsing_bam_filename

get each field from bam filename

=cut
sub parsing_bam_filename {
   my ($self, $irods_bam_file) = @_;

   my ($id_run, $position, $alig_filter, $align_filter2,  $tag_index, $tag_index2, $alig_filter3, $align_filter4 )
                      = $irods_bam_file =~ /\/(\d+)_(\d)(_(\w+))?(\#(\d+))?(_(\w+))?[.]bam$/mxs;

   $align_filter2 ||= $align_filter4;

   return ($id_run, $position, $tag_index2, $align_filter2);
}

=head2 check_bam_header

check_bam_header to see the bam is aligned or not and what reference is used for alignment

=cut
sub check_bam_header {
   my ($self, $irods_bam_file) = @_;

   my $samtools_view_cmd = $self->samtools_irods_cmd . qq[ view -H irods:$irods_bam_file];

   $self->verbose() && $self->log($samtools_view_cmd);

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
      my @ref = grep {/\/(nfs|lustre)\/\S+\/references/mxs} split /\s/mxs, $last_bwa_pg_line;
      if(scalar @ref){
         $reference = $ref[0];
      }
   }

   waitpid $pid, 0;

   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak qq{Failed: $samtools_view_cmd};
   }

   close $sam_header_fh or croak "can not close samtools view output: $ERRNO";

   my @refs = ();
   if($reference){
      push @refs, $reference;
   }

   return {reference => \@refs, alignment => [$sq_tag_found] };
}

=head2 get_irods_file_md5

get md5 value for irods file from ils command

=cut

sub get_irods_file_md5 {
   my ($self, $irods_bam_file) = @_;

   my $irods_loader = npg_common::irods::Loader->new();

   my $md5 = $irods_loader->_get_irods_md5( $irods_bam_file );

   return { md5 => [$md5] };
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

=item npg_common::irods::Loader

=item npg_common::roles::log

=item npg_common::roles::software_location

=item IPC::Open3;

=item POSIX qw(strftime)

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2011 GRL, by Guoying Qi

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
