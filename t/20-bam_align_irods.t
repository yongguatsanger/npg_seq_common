#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-08-02
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Test::More tests => 3;

use_ok('npg_common::bam_align_irods');
{
    my $bam = npg_common::bam_align_irods->new(id_run => 5428, lane => 2, tag_index => 1, working_dir=>'/tmp', task=>'bam_realign');
    is($bam->bam_filename(), '5428_2#1.bam', 'correct bam file name');
    is($bam->old_bam_filename(), '5428_2#1_old.bam', 'correct old bam file name');
}

1;
__END__
