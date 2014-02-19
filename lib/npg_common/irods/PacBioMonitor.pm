#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-02-10
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::PacBioMonitor;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use npg_common::irods::run::PacBio;
use LWP::UserAgent;
use JSON;
use npg_common::irods::Loader;
use Data::Dumper;
use Date::Calc qw(:all);

with qw{MooseX::Getopt npg_common::roles::log};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::PacBioMonitor

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head2 api_uri

pacbio api uri to check Primary Analysis Job Status

=cut
has 'api_uri'           => (isa           => q{Str},
                            is            => q{rw},
                            default       => q{http://pacbio1-1:8081/Jobs/PrimaryAnalysis/Query},
                            documentation => q{PacBio API uri to check Primary Anylysis status},
                           );

=head2 from_date

which date to start

=cut

has 'from_date'         => (isa           => q{Str},
                            is            => q{rw},
                            lazy_build    => 1,
                            documentation => q{which date to start},
                           );
sub _build_from_date {
     my @today = Today();

     my @two_weeks_ago = Add_Delta_Days(@today, -14);## no critic (ProhibitMagicNumbers)

     return sprintf q{%4d-%02d-%02d}, @two_weeks_ago;
}

=head2 before_date

before which date

=cut

has 'before_date'       => (isa           => q{Str},
                            is            => q{rw},
                            documentation => q{before which date},
                           );


=head2 status

which status to check

=cut

has 'status'            => (isa           => q{Str},
                            is            => q{rw},
                            default       => q{Complete},
                            documentation => q{which status to check},
                           );

has 'collection_list'   => (isa           => q{ArrayRef},
                            is            => q{rw},
                            lazy_build    => 1,
                            documentation => q{a list of collections on the instrument},
                           );

sub _build_collection_list {
   my ($self) = @_;

   my $uri = $self->api_uri();
   my  $query = q{status=} . $self->status() . q{&after=} . $self->from_date();
   if( $self->before_date() ){
      $query .= q{&before=} . $self->before_date();
   }
   $uri .= q{?}.$query;

   $self->log("Checking PacBio webservice: $uri");
   my $ua = LWP::UserAgent->new();
   my $response   = $ua->get($uri);

   if ( $response->is_success ) {
     my $list = from_json( $response->content() );
     return $list;
   }

   $self->log("PacBio web service is not available: $uri");

   return [];
}

=head2 process

main method to call to process the collection lsit

=cut

sub process{
  my $self = shift;

  my $collections = $self->collection_list();

  my $count = 0;

  $self->log(scalar  @{$collections}  . q{ runs returned});

  foreach my $collection ( @{$collections} ) {

     $self->log(q{--------------------------});

     my $path = $collection->{OutputFilePath};
     my $collection_order_perwell = $collection->{CollectionOrderPerWell};
     my $collection_number        = $collection->{CollectionNumber};
     my $well                     = $collection->{Well};
     my $plate                    = $collection->{Plate};
     my $index_of_look            = $collection->{IndexOfLook};

     if( !$plate || !$well || !$collection_order_perwell || !$collection_number || !defined $index_of_look){
        $self->log(q{Information not enough for this collection});
        $self->log( qq{\n} . Dumper($collection) );
        next;
     }

     $self->log( "Plate $plate well $well CollectionOrderPerWell $collection_order_perwell CollectionNumber$collection_number IndexOfLook $index_of_look" );

     if( !$path || $path !~ /pacbio1[-]1[.]internal[.]sanger[.]ac[.]uk/mxs ){
         $path ||= q{none};
         $self->log( "IGNORE non-sanger or no path: $path");
         next;
     }else{
        $self->log( $path );
     }

     my $irods = npg_common::irods::Loader->new(file => q{none}, resource=>q{seq});
     my $results = $irods->imeta_query( {type => q{bas},
                           well               => $well,
                           collection_number  => $collection_number,
                           run                => $plate,
                           index_of_look      => $index_of_look,
                          } );
     if(scalar @{$results} ){
         $self->log('This colleciton is already in iRODs');
         next;
     }

     my $local_path = $self->_get_local_path($path);
     if(! -e $local_path ) {
        $self->log("Path not exist: $local_path" );
        next;
     }

     $count++;
     my $collection_folder = $well . q{_} . $collection_order_perwell;
     $local_path =~ s/\/$collection_folder//mxs;

     my $collection_loader = npg_common::irods::run::PacBio
           ->new (
                 runfolder_path    => $local_path,
                 collection_folder => $collection_folder,
                 index_of_look     => $index_of_look,
           );
     $collection_loader->process();
  }

  $self->log("$count collections loaded");

  return 1;
}

sub _get_local_path {
   my ($self, $uri) = @_;

   $uri =~ s/pbids[:]\/\/pacbio1[-]1[.]internal[.]sanger[.]ac[.]uk/\/nfs\/sf45\/pacbio/mxs;

   return $uri;
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

=item npg_common::irods::run::PacBio

=item npg_common::irods::Loader

=item npg_common::roles::log

=item LWP::UserAgent

=item JSON

=item Data::Dumper

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
