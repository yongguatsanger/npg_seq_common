#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2012-05-22
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 13;
use Test::Exception;

use_ok('npg_common::irods::EbiSubmissionData');

{
  my $from_date = '2012-05-19';
  my $stop_date = '2012-05-19';

  my $submission = npg_common::irods::EbiSubmissionData->new(from_date => $from_date,
                                                             stop_date => $stop_date,
                                                             verbose=>1,
                                                            );

  is( $submission->from_date(), $from_date, 'from date set' );
  is( $submission->stop_date(), $stop_date, 'stop date set' );

  throws_ok { $submission->_get_irods_file(1000, 3, undef, "1000_2.bam") } qr/File name 1000_2.bam doesn't match id_run 1000 or lane 3/, 'wrong lane number';
  throws_ok { $submission->_get_irods_file(1000, 3, undef, "2000_3.bam") } qr/File name 2000_3.bam doesn't match id_run 1000 or lane 3/, 'wrong id_run';
  throws_ok { $submission->_get_irods_file(2000, 3, undef, "2000_3#1.bam") } qr/Tag index not given but exists in file name 2000_3#1.bam/, 'Tag index not given but exists in file name ';
  throws_ok { $submission->_get_irods_file(2000, 3, 2, "2000_3#1.bam") } qr/Tag index 2 is wrong in file name 2000_3#1.bam/, 'Tag index is wrong in file name ';
  is($submission->_get_irods_file(2000, 3, 2, "2000_3#2.bam"), '/seq/2000/2000_3#2.bam', 'correct irods full name');
  
  my $current_irods_meta_data = {
       ebi_run_acc => {ERR008027 => 1,},
       ebi_sub_acc => {ERA000191 => 1,}
  };
  my $ebi_submission_data = {
       ebi_run_acc => q{ERR008027},
       ebi_sub_acc => q{ERA000192},
       ebi_sub_date=> q{2012-05-20},
       ebi_sub_md5 => q{dfsdfsdfsd453},
  };
  
  my $imeta_commands = $submission->get_imeta_commands(q{test.bam}, $current_irods_meta_data, $ebi_submission_data);
  is(scalar @{$imeta_commands}, 5, 'correct number of imeta commands');
  is($imeta_commands->[2], q{add -d test.bam ebi_sub_acc "ERA000192"}, 'correct third imeta commands');
  
  my $imeta_commands2 = $submission->get_imeta_commands(q{test.bam}, undef, $ebi_submission_data);
  is(scalar @{$imeta_commands2}, 4, 'correct number of imeta commands when no current meta availabe');
  is($imeta_commands2->[0], q{add -d test.bam ebi_run_acc "ERR008027"}, 'correct first imeta commands');
  is($imeta_commands2->[1], q{add -d test.bam ebi_sub_acc "ERA000192"}, 'correct second imeta commands');

  #$submission->process();

}

1;
__END__
