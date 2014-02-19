#########
# Author:        gq1
# Maintainer:    $Author: gq1 $
# Created:       2010-11-15
# Last Modified: $Date: 2010-12-02 12:30:28 +0000 (Thu, 02 Dec 2010) $
# Id:            $Id: 10-irods-run-BamExtraMeta.t 12048 2010-12-02 12:30:28Z gq1 $
# $HeadURL: svn+ssh://svn.internal.sanger.ac.uk/repos/svn/new-pipeline-dev/useful_modules/branches/prerelease-32.0/t/10-irods-run-BamExtraMeta.t $
#

use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 1;
use Test::Exception;

use_ok('npg_common::irods::run::BamRenameMeta');

{
  my $bam = npg_common::irods::run::BamRenameMeta->new(name=>'sample', old_value=>'161494_A1', new_value=>'UK10K_TW1537988',);

  #my @file_list = @{$bam->file_list()};
  #print "@file_list\n";  
  #$bam->process();
 
  #$bam->process_runlist("t/data/UK10K_TWINS_rename-1.csv");
}

1;
__END__
