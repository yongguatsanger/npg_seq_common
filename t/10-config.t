#########
# Author:        jo3
# Maintainer:    $Author$
# Created:       2010_05_26
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$

use strict;
use warnings;

use File::chdir;
use Test::More tests => 41;
use Test::Deep;
use Test::Exception;
use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/msx; $r; };


my %expected_defaults = (
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
    rdbms     => q{MySQL},
);


use_ok('npg_common::config');

ok( %npg_common::config::DEFAULT, 'Defaults hash is not empty' );

cmp_deeply( \%npg_common::config::DEFAULT, \%expected_defaults,
            'Defaults haven\'t changed since last test' );

local $ENV{HOME} = 't/data/file_system/home';
my $test;

lives_ok { $test = npg_common::config->new() } 'Create object';

throws_ok { $test->dsn() } qr/Searched[ ]for[ ]config[.]ini/msx,
          'Croak if information is lacking';


#
# Basic accessor tests
#


{
    local $ENV{dev} = 'all_present';
    $test = npg_common::config->new( conf_name => 'dsn_config1', );

    is( $test->dbhost(), 'there',       'Use config file dbhost' );
    is( $test->dbport(), '213',         'Use config file dbport' );
    is( $test->dbuser(), 'jeff_bloggs', 'Use config file dbuser' );
    is( $test->dbpass(), 'secret4',     'Use config file dbpass' );
    is( $test->dbname(), 'anything',    'Use config file dbname' );
}


{
    local $ENV{dev} = 'live';
    $test = npg_common::config->new( conf_name => 'dsn_config1', );

    is( $test->dbhost(), $npg_common::config::DEFAULT{dbhost},
        'Fall back to default dbhost' );
    is( $test->dbport(), $npg_common::config::DEFAULT{dbport},
        'Fall back to default dbport' );
    is( $test->dbuser(), $npg_common::config::DEFAULT{dbuser},
        'Fall back to default dbuser' );
    is( $test->dbpass(), $npg_common::config::DEFAULT{dbpass},
        'Fall back to default dbpass' );
    is( $test->dbname(), $npg_common::config::DEFAULT{dbname},
        'Fall back to default dbname' );
}


#
# Tests for _build_file_config
#


{
    lives_ok { $test = npg_common::config->new( conf_name => 'config_1' ) }
             'Create object with file name argument';

    is( $test->dbhost(), 'host_1', 'Find the correct file' );


    $test = npg_common::config->new( conf_name => 'config_2' );
    is( $test->dbhost(), 'host_2', 'Call on Config::Auto to find it' );
}


{
    local $File::chdir::CWD = 't/data/file_system';

    lives_ok {
                $test = npg_common::config->new(
                            conf_name => 'config_1',
                            conf_dir  => 'other_directory',
                )
             }
             'Create object with directory argument';

    is( $test->dbhost(), 'host_1_otherdir', 'Use the argument' );


    $test = npg_common::config->new( conf_name => 'no_such_file_I_hope' );

    throws_ok { $test->dbhost() }
              qr{ Searched[ ]for[ ]\S+\s+
                  Included[ ]paths[ ].+
                  No[ ]config[ ]file[ ]found! }msx,
              'Throw an exception if the config file can\'t be found';

    is(
        ${ $test->warnings() }[-1],
        'Exact match to no_such_file_I_hope not found. Trying variants.',
        'Inspect for warning message'
    );
}


{
    local $ENV{dev} = 'special';

    $test = npg_common::config->new( conf_name => 'config_1' );

    is( $test->dbhost(), 'host_1_envtest', 'Find the correct section' );
}


{
    local $ENV{dev} = 'made_up';
    $test = npg_common::config->new( conf_name => 'config_2', );

    throws_ok { $test->dbhost() } qr/No[ ]made_up[ ]section[ ]found/msx,
              'Croak if section not found';
}


{
    local $ENV{dev} = 'dot_file_only';
    $test = npg_common::config->new( conf_name => 'config_1', );

    $test->dbhost();

    is(
        ${ $test->warnings() }[-1],
        'dot_file_only section not found. Looking in variants.',
        'Inspect for warning message'
    );
}


{
    local $File::chdir::CWD = 't/data/file_system';
    local $ENV{dev}         = 'test';
    $test = npg_common::config->new( conf_name => 'config_1',
                                     conf_dir  => 'other_directory', );

    is( $test->dbhost(), 'host_1_t_test',
        'Modify search path for $' . 'dev eq \'t\'' );
}


#
# Tests for _build_rdbms
#


{
    $test = npg_common::config->new( conf_name => 'rdbms_config1', );
    is( $test->rdbms(), 'no_such_thing', 'Return rdbms from config file' );
}


{
    local $ENV{dev} = 'dev';
    $test = npg_common::config->new( conf_name => 'rdbms_config1', );
    is( $test->rdbms(), $npg_common::config::DEFAULT{rdbms},
        'Otherwise use default if $' . 'dev ne \'test\'' );
}


{
    local $ENV{dev} = 'test';
    $test = npg_common::config->new( conf_name => 'rdbms_config1', );
    is( $test->rdbms(), 'SQLite',
        'Override default if $' . 'dev eq \'test\'' );
}


#
# Tests for _build_dsn
#


{
    $test = npg_common::config->new( conf_name => 'dsn_config1', );
    is( $test->dsn(), 'this_would_never_work',
        'Return a ready-made dsn from the config file' );
}


{
    local $ENV{dev} = 'sqlite_test';
    $test = npg_common::config->new( conf_name => 'dsn_config1',
                                     rdbms     => 'SQLite', );
    like( $test->dsn(), qr/^dbi:SQLite:dbname=\S+/msx, 'Return an SQLite dsn' );
}


{
    local $ENV{dev} = 'no_dbname_test';
    $test = npg_common::config->new( conf_name => 'dsn_config1', );

    throws_ok { $test->dsn() } qr/No[ ]database[ ]defined/msx,
              'Croak on no dbname when asked for dsn';
}


{
    local $ENV{dev} = 'no_dbport_test';
    $test = npg_common::config->new( conf_name => 'dsn_config1', );

    throws_ok { $test->dsn() } qr/No[ ]port[ ]defined/msx,
              'Croak on no dbport when asked for dsn';
}


{
    local $ENV{dev} = 'no_dbhost_test';
    $test = npg_common::config->new( conf_name => 'dsn_config1', );

    throws_ok { $test->dsn() } qr/No[ ]host[ ]defined/msx,
              'Croak on no dbhost when asked for dsn';
}


{
    local $ENV{dev} = 'all_present';
    $test = npg_common::config->new( conf_name => 'dsn_config1',
                                     rdbms     => 'mysql', );

    like( $test->dsn(), qr/^DBI:mysql:database=\S+;host=\S+;port=\d+/msx,
          'Return a MySQL dsn' );


    $test = npg_common::config->new( conf_name => 'dsn_config1',
                                     rdbms     => 'pg', );

    like( $test->dsn(), qr/^dbi:Pg:dbname=\S+/msx, 'Postgres dsn - 1' );


    $test = npg_common::config->new( conf_name => 'dsn_config1',
                                     rdbms     => 'postgres', );

    like( $test->dsn(), qr/^dbi:Pg:dbname=\S+/msx, 'Postgres dsn - 2' );


    $test = npg_common::config->new( conf_name => 'dsn_config1',
                                     rdbms     => 'psql', );

    like( $test->dsn(), qr/^dbi:Pg:dbname=\S+/msx, 'Postgres dsn - 3' );


    $test = npg_common::config->new( conf_name => 'dsn_config1',
                                     rdbms     => 'oracle', );

    like( $test->dsn(), qr/^dbi:Oracle:\S+/msx, 'Return an Oracle dsn' );
}

{
    local $ENV{dev} = 'warehouse3_ro';
    lives_and {
      $test = npg_common::config->new( conf_name => 'config_2',
                                     rdbms     => 'mysql', );
      is $test->dbuser, 'warehouse_ro';
    } 'load a Perl format config';
    cmp_deeply $test->dbattr, {mysql_enable_utf8=>1}, 'read and provide attr hash';
}


#
# Tests for _build_dbh go here
#

1;
