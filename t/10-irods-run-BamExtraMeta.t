#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2010-11-15
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Test::More tests => 2;

use_ok('npg_common::irods::run::BamExtraMeta');
isa_ok( npg_common::irods::run::BamExtraMeta->new(id_run => 5458),
     qw[npg_common::irods::run::BamExtraMeta], 'object test');

1;

