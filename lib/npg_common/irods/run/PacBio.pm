#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-02-08
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::run::PacBio;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use npg_common::irods::Loader;
use XML::LibXML;
use File::Basename;


with qw{MooseX::Getopt npg_common::roles::log};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $DEFAULT_ROOT_DIR   => q{/seq/pacbio};

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::run::PacBio

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

=head1 DESCRIPTION


=head1 SUBROUTINES/METHODS

=head2 runfolder_path 

pacbio runfolder path

=cut
has 'runfolder_path'    => (isa           => q{Str},
                            required      => 1,
                            is            => q{rw},
                            documentation => q{PacBio rundolder path},
                           );

=head2 collection_path 

pacbio collection sub folder name

=cut

has 'collection_folder' => (isa           => q{Str},
                            required      => 0,
                            is            => q{rw},
                            documentation => q{PacBio collection sub folder name},
                           );

=head2 index_of_look 

index of look number for movie

=cut

has 'index_of_look'     => (isa           => q{Int},
                            required      => 0,
                            is            => q{rw},
                            documentation => q{index of look number},
                           );

=head2 file_list

bas and bax h5 file list

=cut

has 'file_list'    => (isa           => q{ArrayRef},
                       is            => q{rw},
                       lazy_build    => 1,
                       documentation => q{bas and bax h5 file list},
                      );

sub _build_file_list {
  my $self = shift;

  my @file_list = ();

  my $collection_subfolder_name = q{*};
  if( $self->collection_folder() ){
     $collection_subfolder_name = $self->collection_folder();
  }

  my $file_pattern = q{*.ba[sx].h5};
  if( defined $self->index_of_look() ){
     $file_pattern = q{*_s} . $self->index_of_look() . q{_} . $file_pattern;
  }
  $file_pattern = $self->runfolder_path().qq{/$collection_subfolder_name/Analysis_Results/$file_pattern};
  $self->log($file_pattern);
  @file_list = glob $file_pattern;

  return \@file_list;
}

sub _get_meta_xml_from_bas_file {
   my ($self, $bas_file_path) = @_;

   $bas_file_path =~ s/Analysis_Results\///mxs;
   $bas_file_path =~ s/bas[.]h5$/metadata\.xml/mxs;
   $bas_file_path =~ s/\b\d+[.]bax[.]h5$/metadata\.xml/mxs; #or PacBioII's extra bax files

   if( ! -e $bas_file_path ){
      croak "The meta xml file does not exist:  $bas_file_path";
   }

   return $bas_file_path;
}

=head2 read_meta_xml

read pacbio metadata xml file and return some of data as a hash ref

=cut

sub read_meta_xml {
   my ($self, $meta_xml) = @_;

   my $irods_meta_data = {};

   my $parser = XML::LibXML->new();
   my $dom = $parser->parse_file( $meta_xml );

   my $run = $dom->getElementsByTagName('Run')->[0];
   my $run_name = $run->getElementsByTagName('Name')->[0]->string_value();
   $irods_meta_data->{run} = $run_name;

   my $sample = $dom->getElementsByTagName('Sample')->[0];
   my $sample_name = $sample->getElementsByTagName('Name')->[0]->string_value();
   $irods_meta_data->{sample} = $sample_name;
   my $well_name = $sample->getElementsByTagName('WellName')->[0]->string_value();
   $irods_meta_data->{well} = $well_name;

   my $instrument_name = $dom->getElementsByTagName('InstrumentName')->[0]->string_value();
   $irods_meta_data->{instrument_name} = $instrument_name;
   my $collection_number = $dom->getElementsByTagName('CollectionNumber')->[0]->string_value();
   $irods_meta_data->{collection_number} = $collection_number;
   my $cell_index = $dom->getElementsByTagName('CellIndex')->[0]->string_value();
   $irods_meta_data->{cell_index} = $cell_index;
   my $set_number = $dom->getElementsByTagName('SetNumber')->[0]->string_value();
   $irods_meta_data->{set_number} = $set_number;

   return $irods_meta_data;
}

sub _get_des_dir {
   my ($self, $local_filename) = @_;

   my ($name,$path,$suffix) = fileparse($local_filename);

   $path =~ s/^.+superfoo/$DEFAULT_ROOT_DIR/mxs;

   return $path;
}

=head2 process

main method to call
get all bas files, load them to irods and add meta data
and load the corresponding metadata xml file as well

=cut

sub process{
  my $self = shift;

  my $bas_file_list = $self->file_list();

  foreach my $bas_file (@{$bas_file_list}){

     my $meta_xml = $self->_get_meta_xml_from_bas_file($bas_file);

     $self->log("---\nLoading $bas_file");
     my $meta_bas = $self->read_meta_xml($meta_xml);
     $meta_bas->{type} = q{bas};
     if( defined $self->index_of_look() ){
        $meta_bas->{index_of_look} = $self->index_of_look();
     }
     my $des_dir_bas = $self->_get_des_dir($bas_file);

     my $loader_bas = npg_common::irods::Loader->new({
       file        => $bas_file,
       collection  => $des_dir_bas,
       meta_data   => $meta_bas,
     });
     $loader_bas->run();


     $self->log("---\nLoading $meta_xml");
     my $meta_meta = {type => q{pacbio meta}};
     my $des_dir_meta = $self->_get_des_dir($meta_xml);
     my $loader_meta = npg_common::irods::Loader->new({
       file        => $meta_xml,
       collection  => $des_dir_meta,
       meta_data   => $meta_meta,
     });
     $loader_meta->run();
  }

  return 1;
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
