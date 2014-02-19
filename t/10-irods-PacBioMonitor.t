#########
# Author:        gq1
# Maintainer:    $Author: gq1 $
# Created:       2011-02-10
# Last Modified: $Date: 2011-02-09 16:00:17 +0000 (Wed, 09 Feb 2011) $
# Id:            $Id: 10-irods-run-PacBio.t 12567 2011-02-09 16:00:17Z gq1 $
# $HeadURL: svn+ssh://svn.internal.sanger.ac.uk/repos/svn/new-pipeline-dev/useful_modules/branches/prerelease-32.0/t/10-irods-run-PacBio.t $
#

use strict;
use warnings;
use Test::More tests => 2;

use_ok('npg_common::irods::PacBioMonitor');

{
  my $monitor = npg_common::irods::PacBioMonitor->new( );
 
  my $pacbio_collection_path = q{pbids://pacbio1-1.internal.sanger.ac.uk/superfoo/Lambda_8_rerun_97/A02_3};
  my $local_path = q{/nfs/sf45/pacbio/superfoo/Lambda_8_rerun_97/A02_3}; 
  is ( $monitor->_get_local_path($pacbio_collection_path), $local_path, 'correct local path');
}

1;
__END__
