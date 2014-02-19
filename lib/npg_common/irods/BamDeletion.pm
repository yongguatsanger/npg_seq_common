#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2012-05-28
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::BamDeletion;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use npg_common::irods::Loader;
use File::Temp qw( tempdir );
use File::Basename;
use File::Spec;

with qw{MooseX::Getopt};
with qw[npg_common::roles::software_location];

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::BamDeletion

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

my $bamDeletion = BamDeletion->new(bam_file => 'irods.bam');

or

my $bamDeletion = BamDeletion->new(bam_file => 'irods.bam', complete_deletion => 1);


$bamDeletion->process();

=head1 DESCRIPTION

Given a bam file, this module will delete all backups of it on irods by default.

If ask for complete deletion, all replications will be deleted but the bam header and meta data will be kept, where meta data type will be changed to bam_header and target will be reset to 0.

If bai file exists for the bam, it will be deleted as well.

=head1 SUBROUTINES/METHODS

=head2 process

main method to call to process

=cut
sub process{
  my $self = shift;

  if( $self->complete_deletion() ){

     #generate a local bam header file from the bam file on irods
     $self->_generate_bam_header();

     #get all the current irods meta data for bam file, and reset the target to 0 and type to bam_header 
     my $bam_header_meta = $self->_get_irods_meta_data();

     #reload bam header with new meta data
     $self->_reload_bam_header($bam_header_meta);

     #remove all replication of bam file    
     $self->_remove_irods_file($self->bam_file());

     #remove all replication of bai file  
     if($self->_bai_file_existed()){
        $self->_remove_irods_file($self->_bai_file());
     }

  }else{

     #remove bam backup
     $self->_remove_backups($self->bam_file());

     #add irods meta data
     $self->_add_backup_removed_meta_data();

     #remove bai backup
     if($self->_bai_file_existed()){
        $self->_remove_backups($self->_bai_file());
     }
  }

  return 1;
}

=head2 working_dir

working directory

=cut
has 'working_dir'              => ( isa           => 'Str',
                                    is            => 'ro',
                                    required      => 0,
                                    default       => sub { tempdir(CLEANUP => 1); },
                                    documentation => 'working directory where temporary files will be created, an optional attribute',
                                  );

=head2 bam_file

irods bam file name

=cut
has 'bam_file'   =>   ( isa           => 'Str',
                        is            => 'ro',
                        required      => 1,
                        documentation => 'irods bam file name',
                      );

=head2 complete_deletion

deletion all copies in irods

=cut
has 'complete_deletion'=> ( isa           => 'Bool',
                            is            => 'ro',
                            required      => 0,
                            documentation => 'boolean option, use it to delete all copies in irods, false by default',
                          );

=head2 rt_ticket

RT ticket number

=cut
has 'rt_ticket'        => ( isa           => 'Int',
                            is            => 'ro',
                            default       => 0,
                            required      => 0,
                            documentation => 'RT ticket number as integer, an optional argument',
                          );

=head2 v

verbose flag

=cut

has 'v'           =>  ( isa           => 'Bool',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'verbose flag',
                      );

has '_header_file_name'        => (  isa           => 'Str',
                                    is            => 'ro',
                                    required      => 0,
                                    lazy_build    => 1,
                                    documentation => 'irods bam header file name',
                                  );
sub _build__header_file_name {
   my $self = shift;
   my $name = basename($self->bam_file);
   $name =~ s/[.]bam$/_header.bam/mxs;
   return $name;
}

=head2 _local_bam_header_file

local bam header file name

=cut
has '_local_bam_header_file'   => ( isa           => 'Str',
                                    is            => 'ro',
                                    required      => 0,
                                    lazy_build    => 1,
                                    documentation => 'local bam header file name',
                                  );
sub _build__local_bam_header_file {
   my $self = shift;
   return  File::Spec->catdir( $self->working_dir(), $self->_header_file_name() );
}

=head2 _irods_bam_header_cmd

command to generate bam header from irods

=cut
has '_irods_bam_header_cmd'   =>  ( isa           => 'Str',
                                    is            => 'ro',
                                    required      => 0,
                                    lazy_build    => 1,
                                    documentation => 'command to generate bam header from irods',
                                  );

sub _build__irods_bam_header_cmd {
   my $self = shift;

   my $bam_file = $self->bam_file();
   my $local_bam_header_file = $self->_local_bam_header_file();
   return $self->samtools_irods_cmd . qq[ view -Hb irods:$bam_file > $local_bam_header_file];
}

has '_bai_file'   =>   ( isa           => 'Str',
                        is            => 'ro',
                        required      => 0,
                        lazy_build    => 1,
                        documentation => 'corresponding bai file name on irods',
                      );
sub _build__bai_file{
   my $self = shift;
   my $bam_file = $self->bam_file();
   $bam_file =~ s/[.]bam$/.bai/mxs;
   return $bam_file;
}

=head2 _remove_backups

remove any obsolete replication and other backups, only keep one copy on the server

If the file and replication numbers not given, use bam's

=cut
sub _remove_backups {
   my ($self, $file) = @_;

   if (!$file) {
      croak 'File should be given';
   }
   my $replication_numbers = $self->_replication_numbers($file);

   $self->_log_message("Remove backups of file $file");

   my $obsolete_replication_nums = $replication_numbers->{obsolete_rep_nums};
   my $replication_nums = $replication_numbers->{rep_nums};

   if( scalar @{$obsolete_replication_nums} ){

      $self->_log_message('There are ' . (scalar @{$obsolete_replication_nums}) . ' obsolete replications to be removed');
      foreach my $rep_num ( @{$obsolete_replication_nums} ) {
         npg_common::irods::Loader->new()->remove_one_replication($file, $rep_num);
      }
   }

   my $num_of_replications =  scalar @{$replication_nums};
   if( $num_of_replications > 1 ){

      $self->_log_message("There are $num_of_replications replications, remove " . ($num_of_replications - 1) . ' backups now' );
      foreach my $i ( 1 .. $num_of_replications -1 ){
         npg_common::irods::Loader->new()->remove_one_replication($file, $replication_nums->[$i]);
      }
   }elsif($num_of_replications == 1){

      $self->_log_message("There is only one copy of file $file, no back up to be removed");
   }else{

      $self->_log_message("There is no copy of file $file");
   }

   return 1;
}

=head2 _add_backup_removed_meta_data

add backup removed meta data

=cut
sub _add_backup_removed_meta_data {
   my ($self) = @_;
   npg_common::irods::Loader->new()->add_meta($self->bam_file(), { backup_removed => 1 });
   return 1;
}

=head2 _generate_bam_header

generate bam header file

=cut
sub _generate_bam_header {
   my $self = shift;

   my $irods_bam_header_cmd = $self->_irods_bam_header_cmd();
   $self->_log_message("Generating bam header file: $irods_bam_header_cmd" );
   if(system $irods_bam_header_cmd){
      croak "Failed to generate bam header: $irods_bam_header_cmd";
   }
   return 1;
}

=head2 _get_irods_meta_data

get all the current irods meta data for bam file

set target = 0 and type = 'bam_header'

=cut
sub _get_irods_meta_data {
   my $self = shift;

   $self->_log_message('Get current bam meta data to generate meta data list for bam header file');
   my $bam_meta_data = npg_common::irods::Loader->new()->_check_meta_data( $self->bam_file() );
   delete $bam_meta_data->{type};
   delete $bam_meta_data->{md5};
   delete $bam_meta_data->{target};

   my %header_meta_data = (target=> 0, type => q{bam_header});
   my $email_meta = 'sample_consent_withdrawn_email_sent';
   if ($self->complete_deletion && $self->rt_ticket && exists $bam_meta_data->{$email_meta}) {
     delete $bam_meta_data->{$email_meta};
     $header_meta_data{$email_meta} = 'RT#' . $self->rt_ticket;
   }

   foreach my $meta_name ( keys %{$bam_meta_data} ){
       push @{$header_meta_data{$meta_name}}, sort keys %{$bam_meta_data->{$meta_name}};
   }
   return \%header_meta_data;
}

sub _reload_bam_header {
   my ($self, $bam_header_meta) = @_;

   $self->_log_message('Loading bam header file to irods');
   my $dir = dirname $self->bam_file;
   return npg_common::irods::Loader->new(
      file       => File::Spec->catfile($self->working_dir, $self->_header_file_name),
      collection => $dir,
      meta_data  => $bam_header_meta,
                                        )->run;
}

sub _remove_irods_file {
   my ($self, $file) = @_;
   $self->_log_message("Remove file $file");
   return npg_common::irods::Loader->new(file => $file)->rm_file();
}

sub _bai_file_existed{
   my $self = shift;
   my $bai_replication_numbers = $self->_replication_numbers($self->_bai_file);
   return @{ $bai_replication_numbers->{rep_nums} } ? 1 : 0;
}

sub _replication_numbers   {
   my ($self, $file) = @_;
   if (!$file) {
     croak 'File should be given';
   }
   return npg_common::irods::Loader->new()->get_replication_numbers($file);
}

sub _log_message {
  my ($self, $message) = @_;
  if ($self->v) {
    warn '***** ' .  $message ."\n";
  }
  return;
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

=item npg_common::roles::software_location

=item File::Temp

=item File::Basename

=item File::Spec

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
