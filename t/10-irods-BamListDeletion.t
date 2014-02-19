#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2012-05-31
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Test::More tests => 6;

use_ok('npg_common::irods::BamListDeletion');

{
  my $deletion = npg_common::irods::BamListDeletion->new();

  isa_ok ($deletion, 'npg_common::irods::BamListDeletion' );
 
  my $expected_iquest_cmd = q{iquest --no-page -z seq "%s/%s" "select COLL_NAME, DATA_NAME where META_DATA_ATTR_NAME = 'ebi_sub_acc' and DATA_NAME not like '%header.bam%'"};
  is( $deletion->_iquest_cmd_for_submitted_bam_files(), $expected_iquest_cmd, 'correct iquest command for submitted files' );

  my $irods_meta = {
          'ebi_sub_acc' => {
                             'ERA135451' => 1
                           },
          'ebi_sub_md5' => {
                             '9c7e74c44380e3798604311042f6a122' => 1
                           },
          'ebi_run_acc' => {
                             'ERR131815' => 1
                           },
          'ebi_sub_date' => {
                              '2012-06-01' => 1
                            },
          'type' => {
                      'bam' => 1
                    },
          'md5' => {
                     '9c7e74c44380e3798604311042f6a122' => 1
                   }
        };
  ok( $deletion->_check_file_deletable('test.bam', $irods_meta), 'bam deletable');
  
  delete $irods_meta->{ebi_sub_acc};
  ok( !$deletion->_check_file_deletable('test.bam', $irods_meta), 'bam not deletable because no ebi_sub_acc');
  
  $irods_meta->{ebi_sub_acc} = {'ERA135451' => 1};
  $irods_meta->{ebi_sub_md5} = {'wrong_values' => 1};
  ok( !$deletion->_check_file_deletable('test.bam', $irods_meta), 'bam not deletable because md5 not match');
}

1;
__END__
