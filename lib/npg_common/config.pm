########
# Author:        jo3
# Maintainer:    $Author$
# Created:       2010-04-28
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$

package npg_common::config;

use Moose;
use MooseX::StrictConstructor;
use Carp;
use Config::Auto;
use DBI;
use English qw(-no_match_vars);
use File::Temp qw(tempfile);
use FindBin;

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };


Readonly::Hash our %DEFAULT => (
    conf_name => q{config.ini},
    conf_dir  => q{data},
    dbhost    => q{},
    dbname    => q{},
    dbpass    => q{},
    dbport    => q{},
    dbuser    => q{},
    dbattr    => undef,
    domain    => q{live},
    dsn       => q{},
    rdbms     => q{MySQL},    # Not the default if domain eq 'test'
);


has [
      'dbhost', 'dbname', 'dbpass', 'dbport',
      'dbuser', 'domain', 'rdbms',  'dsn',
    ] => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
);

has q(dbattr) => (
    is            => 'ro',
    isa           => 'Maybe[HashRef]',
    lazy_build    => 1,
);


# We don't want a lazy build for conf_name because we need to allow
# @arg (tested with accessor) > default
# For conf-Path we want @arg > default
has [ 'conf_name', 'conf_dir' ] => (
    is            => 'ro',
    isa           => 'Str',
);


has 'file_config' => (
    is            => 'ro',
    isa           => 'HashRef',
    lazy_build    => 1,
);


has 'dbh' => (
    is            => 'ro',
    isa           => 'DBI::db',
    lazy_build    => 1,
);


has 'warnings' => (
    is      => 'ro',
    isa     => 'ArrayRef',
    default => sub { return []; },
);


sub _build_dbhost {
    my ($self) = @_;

    return $self->file_config->{dbhost} || $DEFAULT{dbhost};
}


sub _build_dbname {
    my ($self) = @_;

    return $self->file_config->{dbname} || $DEFAULT{dbname};
}


sub _build_dbpass {
    my ($self) = @_;

    return $self->file_config->{dbpass} || $DEFAULT{dbpass};
}


sub _build_dbport {
    my ($self) = @_;

    return $self->file_config->{dbport} || $DEFAULT{dbport};
}


sub _build_dbuser {
    my ($self) = @_;

    return $self->file_config->{dbuser} || $DEFAULT{dbuser};
}

sub _build_dbattr {
    my ($self) = @_;

    return $self->file_config->{dbattr} || $DEFAULT{dbattr};
}


sub _build_domain {
    return ( defined $ENV{dev} ) ? $ENV{dev} : $DEFAULT{domain};
}


sub _build_rdbms {
    my ($self) = @_;

    return $self->file_config->{rdbms} if $self->file_config->{rdbms};

    return 'SQLite' if $self->domain() =~ m/^ test $/imsx;

    return $DEFAULT{rdbms};
}


sub _build_file_config {
    my ($self) = @_;
    my $config;

    my $config_name = $self->conf_name() || $DEFAULT{conf_name};
    my $dir_name    = $self->conf_dir()  || $DEFAULT{conf_dir};
    my $path        = [
                        "./$dir_name",
                        "$ENV{HOME}/$dir_name",
                        "$FindBin::Bin/../$dir_name"
                      ];


    if ( $self->domain() =~ m/^ test $/imsx ) {
        unshift @{ $path }, "./t/$dir_name";
    }


    # First look for the file name 'as is'.
    eval { $config = Config::Auto::parse( $config_name, path => $path );
            1; }

        or do {

            croak $EVAL_ERROR
                if $EVAL_ERROR !~ m/^No[ ]config[ ]file[ ]found!/msx;

            push @{ $self->{warnings} },
                "Exact match to $config_name not found. Trying variants.";

        };

    return $config->{ $self->domain() }
        if defined $config->{ $self->domain() };

    # The above line autovivifies $contig so "$contig" tests true.
    ( scalar keys %{$config} > 0 ) &&
         push @{ $self->{warnings} },
            $self->domain() . ' section not found. Looking in variants.' ;


    # Then try the various Config::Auto expansions.
    local $PROGRAM_NAME = $config_name;
    eval {
        $config = Config::Auto::parse( undef, path => $path );
        1;
    }

    or do {
        croak "Searched for $PROGRAM_NAME"
            . "\nIncluded paths @{$path}\n$EVAL_ERROR";
    };


    croak 'No ' . $self->domain() . ' section found'
        if !defined $config->{ $self->domain() };


    return $config->{ $self->domain() };
}

sub _build_dsn {
    my ($self) = @_;

    return $self->file_config->{dsn} if defined $self->file_config->{dsn};

    my $dsn;
    my $rdbms = lc $self->rdbms();

    if ( $rdbms eq 'sqlite' ) {
        my ( $db_fh, $db_file ) = tempfile( UNLINK => 1 );

        return qq{dbi:SQLite:dbname=$db_file};
    }


    croak 'No database defined' if !$self->dbname();


    if ( $rdbms eq 'mysql' ) {
        croak 'No port defined' if !$self->dbport();
        croak 'No host defined' if !$self->dbhost();

        $dsn = sprintf 'DBI:mysql:database=%s;host=%s;port=%d',
                $self->dbname(), $self->dbhost, $self->dbport();

    }


    if ( $rdbms =~ m/^ postgres | pg | psql $/msx ) {

        $dsn = sprintf 'dbi:Pg:dbname=%s', $self->dbname();

    }


    if ( $rdbms eq 'oracle' ) {

        # Should handle 'sid' here - see DBD::Oracle synopsis.
        $dsn = sprintf 'dbi:Oracle:%s', $self->dbname();

    }

    return $dsn;
}


sub _build_dbh {
    my ($self) = @_;

    my $dbh = DBI->connect( $self->dsn(), $self->dbuser(), $self->dbpass(), $self->dbattr());

    return $dbh;
}


no Moose;

__PACKAGE__->meta->make_immutable;


1;


__END__


=head1 NAME

npg_common::config - flexibly figure out database connection parameters

=head1 VERSION

$Revision$

=head1 SYNOPSIS

    Let the module do all the work:
    C<<my $config = npg_common::config->new();
       my $dbh    = $config->dbh();     # A DBI::db object.
       my $schema = Some::DBIx::Schema->connect( $config->dsn(),
                                                 $config->dbuser(),
                                                 $config->dbpass() );>>

    Give it a hint about what to look for:
    C<<my $config = npg_common::config->new( conf_name => 'proj_config' );>>

    Give it more hints and/or over-ride some parameters:
    C<<my $config = npg_common::config->new( conf_dir => 'etc/settings',
                                             domain   => 'test',
                                             rdbms    => 'SQLite', );>>

    Tell it everything in advance (pointless but illustrative):
    C<<my $config = npg_common::config->new( dbhost => 'mcs123',
                                             dbuser => 'npg_tester',
                                             dbpass => 'secret',
                                             dbport => '3210',
                                             dbname => 'test_db',
                                             rdbms  => 'mysql',
                                             domain => 'live', );>>


    Then get parameters back from it to make a database connection:
    C<<my $username = $config->dbuser();
       my $password = $config->dbpass();
       my $database = $config->dbname();>>

    Or get a dsn string:
    C<<my $dsn = $config->dsn();>>


=head1 DESCRIPTION

This module manages configuration parameters for software projects. At present
it just deals with parameters for database connections, but there is no
obstacle to expanding it for other parameters also.

It defaults to the current practice of the NPD group, expecting to find an
ini-style configuration file called 'config.ini' in a directory called 'data',
together with the use of the environment variable '$dev' (referred to as the
'domain' in this documentation) to choose sections within the configuration
file. If $dev is unset, a default value of 'live' is assumed. The module uses
the functionality of L<Config::Auto|Config::Auto> (q.v.) to make best guesses
about the name, location and format of the configuration file if this
practice does not apply.

It is not absolutely dependant on finding the file and will accept some or all
necessary parameters as constructor arguments to fill in for missing
parameters or to over-ride existing ones. However it will croak if it
discovers that it cannot find any value for a parameter it needs.

When an npg_common::config object has been instantiated, it can be used to
return individual configuration parameters, a dsn string or a database
connection in the form of a DBI::db object.

For some parameters a default value is hard-coded into the module. These
generally fall back to a temporary test database, but can be inspected via
%ngp_common::config::DEFAULT. They are over-ridden if a configuration file is
found. Those configuration file parameters, in turn, are over-ridden by
arguments passed to the constructor.

=head2 Constructor Arguments

=over

=item B<conf_name>

The base name of the configuration file. See the documentation for
L<Config::Auto|Config::Auto> for information on how it builds on this base.

If none is supplied the module falls back on the default value. This default
is set to reflect current practice in our group and the module will search
only for an exact match to the value in this case.

=item B<conf_dir>

A directory name to use when adding to the list of paths in which to search
for the config file. The paths added will be ./$conf_dir, ~/$conf_dir and
$FindBin::Bin/../$conf_dir and they are added after the current directory, but
before Config::Auto's list).

If the domain is 'test' the path 't/$conf_dir' is unshift()ed to the front of
the above list.

=item B<domain>

The configuration domain name. This is usually 'live', 'dev', or 'test', but
can be anything that matches a section of the configuration file. Tests on the
value of this attribute are case-insensitive.

=item B<rdbms>

The database management system. The currently supported, and case-independent,
values are 'sqlite', 'mysql', 'psql', (also 'pg' and 'postgres'), and
'oracle'.

=item B<dbname>

The database name.

=item B<dbuser>

The database username.

=item B<dbpass>

The database password.

=item B<dbhost>

The database host.

=item B<dbport>

The database port.

=item B<dsn>

The dsn string consisting of some or all of the above elements and formed
according to the value of the 'rdbms' parameter.

=back

=head2 Accessors

Named as the constructor arguments above with some additional ones.

=over

=item B<schema>

A DBIx::Schema database connection object.

=item B<dbh>

A DBI::db database connection handle.

=back

=head1 SUBROUTINES/METHODS

None.

=head1 CONFIGURATION AND ENVIRONMENT

Config files are sought and parsed using L<Config::Auto|Config::Auto> so the format and location of the file is very flexible. See the documentation of
that module for the detail of how it operates.

The 'domain' parameter referred to by this module corresponds to the names of
sections of the config file. If not passed to the module constructor, the
environment variable '$dev' is used to set this attribute. If $dev is not set
a default is used.

=head1 INCOMPATIBILITIES

None known.

=head1 DIAGNOSTICS

=head1 DEPENDENCIES

=over

=item Moose

=item DBI

=item Readonly

=item Carp

=item English

=item Config::Auto

=item FindBin

=item npg_tracking::Schema

=item File::Temp

=back

=head1 BUGS AND LIMITATIONS

No methods for DBI's \%attr parameter (AutoCommit, etc).

The module performs two searches. One for an exact match to the config file
name and then, if no file is found, another for variants on the file name as
per Config::Auto's documentation. The same list of directories is searched
each time.

A matching config file can be found in a number of places in that list. The
module will read the first one it finds, and will then return. In a situation
where a matching file is found but the particular section required is not
present, the module will not look any further along its search path. (However
it will try the second search path.)

This can give rise to a situation where the module can return a configuration
that it obtained from a file other than the one the user expected.

The user should check for this possibility by inspecting the arrayref returned
by $config->warnings(). If it is empty, the module found the required section
in the first matching config file that it found.

If the required section is not found in any file, the module croaks.


Please contact the author with any other bugs/limitations found.

=head1 AUTHOR

John O'Brien, E<lt>jo3@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by John O'Brien

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

