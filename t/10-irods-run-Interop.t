#########
# Author:        Jennifer Liddle (js10)
# Maintainer:    $Author$
# Created:       2012-02-07
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Test::More tests => 3;

use_ok('npg_common::irods::run::Interop');

my $EXIST_EXECUTABLES = `which icd`;

# ensure we're in the user's iRODs home
my $have_irods;
if ($EXIST_EXECUTABLES) {
  `icd`;
  $have_irods = system("imkdir test") == 0 ? 1 : 0;
}

SKIP: {
  unless ($EXIST_EXECUTABLES) {
    skip 'unable to access iRODS executables', 2;
  }
  elsif (!$have_irods) {
    skip 'unable to create iRODS test dir (try kinit to log in to iRODS)', 2;
  }

  my $interop = new npg_common::irods::run::Interop(id_run => 5174, 
	      runfolder_path=>'t/data/staging/120118_HS22_07385_A_D0GKHACXX',
	      irods_root=>'test/test/',
	                                                 );
  is(ref $interop, 'npg_common::irods::run::Interop');
  is($interop->process(),1);

  `irm -r test`;
}

1;

