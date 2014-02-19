#############
# $Id$
# Created By: ajb
# Last Maintained By: $Author$
# Created On: 2010-09-22
# Last Changed On: $Date$
# $HeadURL$

package npg_common::roles::run::intensities::config;
use Moose::Role;
use XML::LibXML;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$LastChangedRevision$ =~ /(\d+)/mxs; $r; };
## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::roles::run::intensities::config

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

  package My::Package;
  use Moose;

  ...something which provides method intensity_path...

  with qw{npg_common::roles::run::intensities::config};

=head1 DESCRIPTION

This role reads information out of the config.xml file in the intensities directory.
You probably want to have loaded npg_tracking::illumina::run::folder first, in order
to have the intensity_path method, although you can provide your own.

The method intensities_config_xml provides you with an XML::LibXML::Document object of
the file, which you can use as you see fit. Other methods are helper methods to get the
information.

=head1 SUBROUTINES/METHODS

head2 intensities_config_xml

Provides an XML::LibXML::Document object of the config.xml document from the Intensities directory.
Other methods are helper methods to obtain the information 

  my $oConfigXml = $oMyPackage->intensities_config_xml();

=cut

has q{intensities_config_xml} => (
  isa => q{Object},
  is  => q{ro},
  init_arg => undef,
  lazy_build => 1,
);

sub _build_intensities_config_xml {
  my ( $self ) = @_;

  return XML::LibXML->load_xml(
    location => $self->intensity_path() . q{/config.xml},
  );
}

=head2 first_cycle

The declared first cycle number

  my $iFirstCycle = $oMyPackage->first_cycle();

=head2 last_cycle

The declared last cycle number

  my $iLastCycle = $oMyPackage->last_cycle();

=cut

has q{first_cycle} => (
  isa => q{Int},
  is  => q{ro},
  lazy_build => 1,
  writer => q{_set_first_cycle},
);

sub _build_first_cycle {
  my ( $self ) = @_;
  return $self->_cycles( q{first_cycle} );
}

has q{last_cycle} => (
  isa => q{Int},
  is  => q{ro},
  lazy_build => 1,
  writer => q{_set_last_cycle},
);

sub _build_last_cycle {
  my ( $self ) = @_;
  return $self->_cycles( q{last_cycle} );
}

has q{rta_version} => (
  isa => q{Str},
  is  => q{ro},
  lazy_build => 1,
  init_arg => undef,
  writer => q{_set_rta_version},
);

sub _build_rta_version {
  my ( $self ) = @_;
  return $self->_rta_info( q{rta_version} );
}

has q{rta_name} => (
  isa => q{Str},
  is  => q{ro},
  lazy_build => 1,
  init_arg => undef,
  writer => q{_set_rta_name},
);

sub _build_rta_name {
  my ( $self ) = @_;
  return $self->_rta_info( q{rta_name} );
}

#############
# private methods

# parses out all cycle numbers, and populates them into all, so that a build of one, builds all
sub _cycles {
  my ( $self, $return_key ) = @_;

  my $run_el      = ( $self->intensities_config_xml()->getElementsByTagName( 'Run' ) )[0];
  my $cycles_el   = ($run_el->getChildrenByTagName( 'Cycles' ))[0];
  my $cycle_first = $cycles_el->getAttribute( 'First');
  my $cycle_last  = $cycles_el->getAttribute( 'Last' );

  $self->_set_first_cycle( $cycle_first );
  $self->_set_last_cycle( $cycle_last );

  my $cycles = {
    first_cycle => $cycle_first,
    last_cycle  => $cycle_last,
  };

  return $cycles->{ $return_key };
}

# parses out the rta version and name
sub _rta_info {
  my ( $self, $return_key ) = @_;

  my $ele = $self->intensities_config_xml()->getElementsByTagName('Software')->[0];
  my $version = $ele->getAttribute('Version');
  my $name = $ele->getAttribute('Name');

  $self->_set_rta_version( $version );
  $self->_set_rta_name( $name );

  my $rta_info = {
    rta_version => $version,
    rta_name    => $name,
  };

  return $rta_info->{ $return_key };
}

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item Carp

=item English -no_match_vars

=item Readonly

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

$Author$

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by Andy Brown (ajb@sanger.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
