#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-07-01
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 9;
use Test::Exception;
use JSON;
use File::Temp qw(tempdir);
use File::Which;

use_ok('npg_common::irods::BamMetaUpdater');

# create a place holder if samtools_irods not visible on the path 
my $samtools_irods = which(q[samtools_irods]);
my $samtools_irods_is_executable = 1;
unless ($samtools_irods && -s $samtools_irods) {
   my $working_dir = tempdir(CLEANUP => 1);
   $samtools_irods = qq[$working_dir/samtools_irods];
   `touch $samtools_irods`;
   `chmod 755 $samtools_irods`;
} 

{
  my $bam = npg_common::irods::BamMetaUpdater->new(
                                   verbose => 1,
                                   id_run => 6345,
                                   lane => 6,
                                   tag_index => 6,
                                   samtools_irods_cmd => $samtools_irods,
                                );

  ok($bam->samtools_irods_cmd, 'has samtools_irods command');

  
  my $result = [$bam->parsing_bam_filename("/seq/6345/6345_5.bam")];
  is_deeply( $result, [6345, 5, undef, undef], 'id_run, lane' );

  $result = [$bam->parsing_bam_filename("/seq/6345/6345_5_phix.bam")];
  is_deeply( $result, [6345, 5, undef, 'phix'], 'id_run, lane, phix' );

  $result = [$bam->parsing_bam_filename("/seq/6345/6345_5_phix#6.bam")];
  is_deeply( $result, [6345, 5, 6, 'phix'], 'id_run, lane, phix, plex' );

  $result = [$bam->parsing_bam_filename("/seq/6345/6345_5_nonhuman#6.bam")];
  is_deeply( $result, [6345, 5, 6, 'nonhuman'], 'id_run, lane, nohuman plex' );
  
  $result = [$bam->parsing_bam_filename("/seq/6345/6345_5#6.bam")];
  is_deeply( $result, [6345, 5, 6, undef], 'id_run, lane, plex' );
  
  $result = [$bam->parsing_bam_filename("/seq/6345/6345_5#6_phix.bam")];
  is_deeply( $result, [6345, 5, 6, 'phix'], 'id_run, lane, plex pf0' );
  
  is($bam->_remove_history_date("[2012-06-20T15:38:01] 3836??3800T3"), "3836??3800T3", 'remove history date');
}

1;
__END__
