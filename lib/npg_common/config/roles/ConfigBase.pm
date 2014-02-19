#############
# $Id$
# Created By: ajb
# Last Maintained By: $Author$
# Created On: 2010-09-08
# Last Changed On: $Date$
# $HeadURL$

package npg_common::config::roles::ConfigBase;
use Moose::Role;
use Config::Auto;
use Config::Any;
use English qw(-no_match_vars);
use File::Temp qw(tempfile);
use FindBin;
use Carp;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$LastChangedRevision$ =~ /(\d+)/mxs; $r; };

## no critic (Documentation::RequirePodAtEnd)

Readonly::Hash our %DEFAULT => (
    CONF_NAME => q{config.ini},
    CONF_DIR  => q{data},
    DOMAIN    => q{live},
);

=head1 NAME

npg_common::config::roles::ConfigBase

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

  package My::Config::Class;
  use Moose;
  with qw{npg_common::config::roles::ConfigBase};
  ....
  no Moose;
  1;

=head1 DESCRIPTION

This role sets up the base for locating and reading in a Config file, utilising Config::Auto, so that it should
be able to use any of a range of config formats. (See the cpan page for Config::Auto for this)

=head1 SUBROUTINES/METHODS

=head2 conf_name

accessor for the conf_name. can be set on construction or else defaults to config.ini

=head2 conf_dir

accessor for the conf_dir. can be set on construction or else defaults to data

=head2 domain

accessor for the environment domain to retrieve from the config file
will build itself from $ENV{TEST} or the default value (in that preferential order)

default = live

=cut

has q{conf_name} => (
  is            => q{ro},
  isa           => q{Str},
  documentation => q{Configuration filename. Default: } . $DEFAULT{CONF_NAME},
  default       => $DEFAULT{CONF_NAME},
);

has q{conf_dir} => (
  is            => q{ro},
  isa           => q{Str},
  documentation => q{Configuration directory. Default: } . $DEFAULT{CONF_DIR},
  default       => $DEFAULT{CONF_DIR},
);

has q{domain} => (
  is            => q{ro},
  isa           => q{Str},
  documentation => q{Configuration filename. Default: } . $DEFAULT{DOMAIN},
  lazy_build    => 1,
);

=head2 conf_path

accessor for the conf_path. can be set on construction or else defaults to current working directory
this path will be used as is, it will not add conf_dir to the end

=cut

has q{conf_path} => (
  is            => q{ro},
  isa           => q{Str},
  documentation => q{path to directory containing configuration files. Default: Current working Directory},
  default       => q{.},
);

=head2 croak_config_read_if_warnings

Boolean flag to tell the config reader to croak if there are any warnings whilst
trying to read in the config file

This is deliberately rw, so that it can be set/changed dynamically be class internal methods

=cut

has q{croak_config_read_if_warnings} => (
  isa => q{Bool},
  is  => q{rw},
  documentation => q{Boolean. Default: 0},
);

=head2 file_config

accessor contains the hash of details taken from the config file that has been read
use this to obtain the values which you should then 'build' into your accessors

=cut

has q{file_config} => (
    is            => q{ro},
    isa           => q{HashRef},
    lazy_build    => 1,
);

has q{config_warnings} => (
  traits        => ['Array'],
  isa           => q{ArrayRef[Str]},
  is            => q{ro},
  init_arg      => undef,
  default       => sub { [] },
  handles       => {
    all_config_warnings   => q{elements},
    no_config_warnings    => q{is_empty},
    count_config_warnings => q{count},
    first_config_warning  => q{shift},
    add_config_warning    => q{push},
  },
);

################
# private methods

sub _build_domain {
  my ( $self ) = @_;
  return $ENV{TEST} || $DEFAULT{DOMAIN};
}

# reads the config file and utilising domain, creates the file_config hashref
sub _build_file_config { ##no critic (Subroutines/ProhibitExcessComplexity)
  my ( $self ) = @_;

  my $config;

  my $config_name = $self->conf_name();
  my $dir_name    = $self->conf_dir();
  my $path        = [
                      $self->conf_path(),
                      "./$dir_name",
                      "$ENV{HOME}/$dir_name",
                      "$FindBin::Bin/../$dir_name"
                    ];


  if ( $self->domain() =~ m/^ test $/imsx ) {
      unshift @{ $path }, "./t/$dir_name";
  }

  # First look for the file name 'as is'.    
  eval {
#diag (Config::Any->extensions());
#    $config = Config::Auto::parse( $config_name, path => $path );
    foreach my $path_dir ( @{ $path } ) {

      my $file_name_and_path = $path_dir . q{/} . $config_name;

      $config = Config::Any->load_files( { files => [ $file_name_and_path ], use_ext => 1, } );

      if ( scalar @{ $config } ) {
        $config = $config->[0]->{ $file_name_and_path };
        last;
      }
    };

    if ( ref $config eq q{ARRAY} && ! scalar @{ $config } ) {
      croak q{No config file found!};
    };

    $config;
  } or do {

    $config = undef; # we want rid of the arrayref;
    # if the error is only that it didn't find the config file, don't croak just yet
    if ( $EVAL_ERROR =~ m/\ANo[ ]config[ ]file[ ]found!/msx ) {
      $self->add_config_warning( qq{Exact match to $config_name not found. Trying variants.} );
    } else {
      croak $EVAL_ERROR;
    }

  };

  # if $config is not set, don't accidently autovivify
  # return now if you have what you are looking for
  if ( defined $config and defined $config->{ $self->domain() } ) {
    return $config->{ $self->domain() };
  }

  # if a config was obtained, then add a warning to say that it was found, but no domain located
  $config && $self->add_config_warning( $self->domain() . q{ section not found. Looking in variants.} );

  # string together warnings so that they can be printed with any future croaks
  my $warnings_string = q{};
  if ( $self->count_config_warnings() ) {

    $warnings_string = join qq{\n\t}, $self->all_config_warnings();
  }

  # if the user has specified croak if warnings are obtained, then do just that
  # most likely if you don't want it trying to find any other possible config file
  if ( $self->croak_config_read_if_warnings() &&  $self->count_config_warnings() ) {
    croak qq{Config warnings found - croaking:\n\t$warnings_string\n};
  }

  # Then try the various Config::Auto expansions.
  local $PROGRAM_NAME = $config_name;
  eval {
    $config = Config::Auto::parse( undef, path => $path );

  } or do {

    # now just croak out regardless
    my $msg = q{};
    if ( $self->count_config_warnings() ) {
      $msg = qq{WARNINGS PRIOR TO CURRENT CROAK:\n\t$warnings_string\n};
    }

    croak $msg . qq{CROAK:\n\tSearched for $PROGRAM_NAME\nIncluded paths @{$path}\n$EVAL_ERROR};

  };

  # if the domain cannot be found in the config, then croak
  if ( ! defined $config->{ $self->domain() } ) {
    my $msg = q{};
    if ( $self->count_config_warnings() ) {
      $msg = qq{WARNINGS PRIOR TO CURRENT CROAK:\n\t$warnings_string\n};
    }
    croak $msg . q{No } . $self->domain() . q{ section found};
  }

  # if the class has a log method, and there are config warnings, then log them for reference
  if ( $self->can( q{log} ) && $self->count_config_warnings() ) {
    $self->log( qq{WARNINGS: $warnings_string} );
  }

  return $config->{ $self->domain() };
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

=item Config::Auto

=item File::Temp

=item FindBin

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
