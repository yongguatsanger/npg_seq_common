#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2012-05-31
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::BamListDeletion;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use npg_common::irods::Loader;
use npg_common::irods::BamDeletion;
use IPC::Open3;
use MIME::Lite;

with qw{MooseX::Getopt npg_common::roles::log};
with qw{npg_common::roles::software_location};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $EXIT_CODE_SHIFT => 8;
Readonly::Scalar our $RT_TICKET_EMAIL_ADDRESS => q{new-seq-pipe@sanger.ac.uk};

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::BamListDeletion

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

my $bamListDeletion = BamListDeletion->new();

or

my $bamListDeletion = BamListDeletion->new( bam_list_file => $bam_filename_list_file );

or

my $bamListDeletion = BamListDeletion->new( consent_withdrawn_deletion => 1 );


$bamListDeletion->process();

In the first two cases you can always pass complete_deletion in, otherwise it only remove any backups 

=head1 DESCRIPTION

This module will delete bam files in irods based on a list of bam files in a given file, where only the file names are listed without the irods collection names. For example, 1234_5#6.bam, 1234 will be used as id_run and then the whole bam file path will be /seq/1234/1234_5#6.bam

If the list file is not given, it will query irods to get a list file automatically.

There are two levels of deletions: only remove backups if complete_deletion not given.

The module also can identify bam files with sample_consent_withdraw meta data and send an email to new-seq-pipe@sanger.ac.uk and mark them with sample_consent_withdraw_email_sent.

TODO: set sample_consent_withdraw meta data in irods_meta_data updater when switching to warehouse 3.

=head1 SUBROUTINES/METHODS

=head2 process

main method to call to process

=cut
sub process{
  my $self = shift;

  if($self->bam_list_file()){
      $self->process_bam_list_from_file();
  }elsif( !$self->consent_withdrawn_deletion() ){
      $self->process_submitted_files();
  }

  return 1;
}

=head2 complete_deletion

deletion all copies in irods

=cut
has 'complete_deletion'=> ( isa           => 'Bool',
                            is            => 'ro',
                            required      => 0,
                            documentation => 'deletion all copies in irods',
                          );

=head2 bam_list_file

irods bam list file name

=cut
has 'bam_list_file'   =>   ( isa           => 'Str',
                             is            => 'ro',
                             required      => 0,
                             documentation => 'irods bam list file name',
                           );

=head2 _iquest_cmd_for_submitted_bam_files

iquest to get the submitted irods bam file list

=cut
has '_iquest_cmd_for_submitted_bam_files'  => (
                            isa       => 'Str',
                            is        => 'ro',
                            required  => 0,
                            default   => q{iquest --no-page -z seq "%s/%s" "select COLL_NAME, DATA_NAME where META_DATA_ATTR_NAME = 'ebi_sub_acc' and DATA_NAME not like '%header.bam%'"},
                              );

=head2 process_submitted_files

process ebi submitted files

=cut
sub process_submitted_files {
  my $self = shift;

  my $iquest_cmd = $self->_iquest_cmd_for_submitted_bam_files();

  $self->log($iquest_cmd);

  my $pid = open3( undef, my $iquest_out_fh, undef, $iquest_cmd);

  my $count = 0;
  my $count_deletable = 0;

  while (my $line = <$iquest_out_fh> ) {

     chomp $line;

     if($line =~ /^ERROR/mxs){
         $self->log($line);
         last;
     }elsif($line =~/[.]bam$/mxs){

         $self->log($line);
         $count++;

         if( $self->_check_file_deletable($line) ){

            $count_deletable++;
            my $bam_deletion = npg_common::irods::BamDeletion->new(
                                           bam_file          => $line,
                                           complete_deletion => $self->complete_deletion(),
                                           $self->resolved_paths,
                                                                 );
            $bam_deletion->process();
         }else{
            $self->log("$line is not deletable");
         }
     }
  }

  waitpid $pid, 0;

  if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     carp "Failed: $iquest_cmd";
  }

  close $iquest_out_fh or croak "can not close iquest command output: $ERRNO";

  $self->log("$count_deletable out of $count bam files are deleted");

  return 1;
}

=head2 process_bam_list_from_file

process bam list from a file

=cut
sub process_bam_list_from_file{
   my ($self, $list_file) = @_;

   if(!$list_file){
      $list_file = $self->bam_list_file();
   }

   open my $list_fh, '<', $list_file or croak "can not open file $list_file";## no critic (InputOutput::RequireBriefOpen);

   my $count = 0;
   my $count_deletable = 0;

   while(my $line = <$list_fh> ){

      chomp $line;
      $self->log($line);
      $count++;

      my ($id_run) = $line =~ /(^\d+)_.+[.]bam$/mxs;
      if($id_run){

          my $irods_file = "/seq/$id_run/$line";
          if( $self->_check_file_deletable($irods_file) ){

            $count_deletable++;
            my $bam_deletion = npg_common::irods::BamDeletion->new(
                                           bam_file          => $irods_file,
                                           complete_deletion => $self->complete_deletion(),
                                           $self->resolved_paths,
                                         );
            $bam_deletion->process();
         }else{
            $self->log("$irods_file is not deletable");
         }
      }else{
          $self->log("no id run found from file name $line");
      }
   }

   close $list_fh or croak "can not close file $list_file";

   $self->log("$count_deletable out of $count bam files are deleted");

   return 1;
}

sub _check_file_deletable {
   my ($self, $irods_file, $irods_meta) = @_;

   if(! $irods_meta){
      $irods_meta = npg_common::irods::Loader->new()->_check_meta_data($irods_file);
   }

   my $deletable = 1;

   my @existed_meta_data = qw(ebi_sub_acc ebi_run_acc ebi_sub_md5 ebi_sub_date);
   foreach my $meta (@existed_meta_data){
      if( !$irods_meta->{$meta}
         || ( $irods_meta->{$meta} && scalar keys %{$irods_meta->{$meta}} == 0 )
        ){

         $self->log("$irods_file: $meta meta data does not exist in irods");
         $deletable = 0;
      }
   }

   my $md5;
   my $ebi_sub_md5;

   if($irods_meta->{md5} && scalar keys %{$irods_meta->{md5}} == 1 ){
      ($md5) = keys %{$irods_meta->{md5}}
   }
   if( $irods_meta->{ebi_sub_md5} && scalar keys %{$irods_meta->{ebi_sub_md5}} == 1 ){
      ($ebi_sub_md5) = keys %{$irods_meta->{ebi_sub_md5}}
   }

   if( !$md5 || !$ebi_sub_md5 || $md5 ne $ebi_sub_md5 ){
       $self->log("$irods_file: EBI submission md5 and current irods md5 are different");
       $deletable = 0;
   }

   return $deletable;
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

=item npg_common::irods::BamDeletion

=item npg_common::roles::log

=item npg_common::roles::software_location

=item IPC::Open3

=item MIME::Lite

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2012 GRL, by Guoying Qi

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
