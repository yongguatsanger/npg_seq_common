# $Id$

use strict;
use warnings;
use Test::More tests => 8;
use Test::Exception;
use lib qw{t/lib};
use npg::api::util;
use npg_testing::intweb qw(npg_is_accessible);

local $ENV{'NPG_WEBSERVICE_CACHE_DIR'} = q[t/data/sequence];
local $ENV{dev} = q{dev}; # ensure that we always have the dev server for this test

BEGIN {
  use_ok( q{npg_common::roles::run::status} );
}

use npg_common::role_tests::run_status_no_extra_methods;
use npg_common::role_tests::run_status_has_run_method;

{
  my $status;
  lives_ok {
    $status = npg_common::role_tests::run_status_no_extra_methods->new();
  } q{obtain object ok};

  throws_ok {
    $status->update_run_status( {} );
  } qr{no[ ]id_run[ ]and/or[ ]status_desc[ ]provided}, q{no id_run};

  throws_ok {
    $status->update_run_status( { id_run => 1234, } );
  } qr{no[ ]id_run[ ]and/or[ ]status_desc[ ]provided}, q{no status_desc};
}

SKIP: {
  skip q{Unable to contact dev webserver}, 4 unless npg_is_accessible($npg::api::util::DEV_BASE_URI); 
 {
  local $ENV{'NPG_WEBSERVICE_CACHE_DIR'} = undef;
  my $status;
  lives_ok {
    $status = npg_common::role_tests::run_status_no_extra_methods->new();
  } q{obtain object ok};

  lives_ok {
    $status->update_run_status( {
      id_run => 1234,
      status_desc => q{run pending},
    } );
  } q{status updated};
 }
 {
  my $status;
  lives_ok {
    $status = npg_common::role_tests::run_status_has_run_method->new();} q{obtain object ok};

  lives_ok {
    $status->update_run_status( {
      status_desc => q{run pending},
    } );
  } q{status updated};
 }
}
1;
