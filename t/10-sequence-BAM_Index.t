#########
# Author:        gq1
# Created:       2010-10-05
#

use strict;
use warnings;
use Test::More tests => 3;
use Cwd qw(abs_path);

use_ok('npg_common::sequence::BAM_Index');
{
  my $jar = 't/bin/aligners/picard/picard-tools-1.31/MarkDuplicates.jar';
  my $bam = npg_common::sequence::BAM_Index->new(
                 input_bam     => '/tmp/5233/5233_1#0.bam',
                 bam_index_jar_file => $jar,
               );
  my $java_cmd = $bam->java_cmd;
  is( $bam->output_bam_index(), '/tmp/5233/5233_1#0.bai', 'output bam index file name');
  my $abs_jar_path = abs_path($jar);
  is( $bam->bam_index_cmd(), qq{$java_cmd -Xmx2000m -jar $abs_jar_path INPUT=/tmp/5233/5233_1#0.bam VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'}, 'bam indexing cmd with absolute Picard path');
}

1;
