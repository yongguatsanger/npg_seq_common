#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2012-05-28
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Test::More tests => 31;
use Test::Exception;
use File::Temp qw( tempdir );
use Cwd;

use_ok('npg_common::irods::BamDeletion');
isa_ok( npg_common::irods::BamDeletion->new(
                                   bam_file => '/seq/1000/1000_5#4_test.bam',
                           ), 'npg_common::irods::BamDeletion');

my $working_dir = tempdir(CLEANUP => 1);
my @path = split '/', $working_dir;

my $EXIST_EXECUTABLES = exist_irods_executables();
my $IRODS_TEST_AREA = '/seq/npg/test/' . pop @path;
my $irods_test_area_created = $EXIST_EXECUTABLES ? create_irods_test_area() : 0;

my $expected_bam_header_meta = {
          'target' => {'0' => 1},
          'type'   => {'bam_header' => 1},
          'sample' => {'s1' => 1,'s2' => 1,},
          'md5' => {'c2bb1d4d0523e29714df101f2e410a31' => 1},
          'id_run' => {'1000' => 1},
};

my $samtools_irods = qq[$working_dir/samtools_irods];
`touch $samtools_irods`;
`chmod 755 $samtools_irods`;
 

{
  my $deletion = npg_common::irods::BamDeletion->new(
                          bam_file    => '/seq/1000/1000_5#4_test.bam',
                          working_dir => $working_dir,
                          samtools_irods_cmd => $samtools_irods,
                        );
  is( $deletion->_header_file_name(), '1000_5#4_test_header.bam', 'correct irods bam header file name' );
  is( $deletion->_local_bam_header_file(), $working_dir . q{/1000_5#4_test_header.bam}, 'correct local bam header file name' );

  is( $deletion->_irods_bam_header_cmd(), "$samtools_irods view -Hb irods:/seq/1000/1000_5#4_test.bam > $working_dir/1000_5#4_test_header.bam", 'correct irods bam header command');

  is( $deletion->_bai_file(), '/seq/1000/1000_5#4_test.bai', 'correct bai file name');
}


SKIP: {

  unless ($EXIST_EXECUTABLES) {
    skip 'unable to access iRODS executables', 25;
  }
  elsif (! $irods_test_area_created) {
    skip 'unable to create iRODS test area (try kinit to log in to iRODS)', 25;   }

{
  my $loader = new npg_common::irods::Loader(
             file       => 't/data/test.bam',
             collection => $IRODS_TEST_AREA,
             meta_data  => {type      => 'bam',
                            id_run    => '1000',
                            sample    => [qw(s1 s2)],
                            target    => 1,
                           },
  );
  $loader->run();

  my $bam = $IRODS_TEST_AREA . '/test.bam';

  my $deletion;

  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run', 1 unless ($ENV{TOOLS_INSTALLED});
    $deletion = npg_common::irods::BamDeletion->new(
                 bam_file    => $bam,
                 working_dir => $working_dir,
                 v => 1,
                 #samtools_irods_cmd => $samtools_irods,
                );
    $deletion->samtools_irods_cmd();
  } 

  SKIP: {
    skip 'samtools_irods executable found', 1 if ($ENV{TOOLS_INSTALLED});
    $deletion = npg_common::irods::BamDeletion->new(
                 bam_file    => $bam,
                 working_dir => $working_dir,
                 v => 1,
                 samtools_irods_cmd => $samtools_irods,
                );
  }

  my $header = $IRODS_TEST_AREA . '/' . $deletion->_header_file_name;

  is_deeply($deletion->_replication_numbers($deletion->bam_file), {rep_nums => [0, 1], obsolete_rep_nums => [] },
      'correct replication numbers');

  is_deeply($deletion->_replication_numbers($deletion->_bai_file),
      {rep_nums => [], obsolete_rep_nums => [] }, 'correct replication numbers for bai file');

  ok(!$deletion->_bai_file_existed(), 'no bai file');

  lives_ok {$deletion->_generate_bam_header()} 'generate bam header';

  my $expected_meta  = {
          'target' => 0,
          'sample' => ['s1', 's2'],
          'type' => 'bam_header',
          'id_run' => ['1000'],
        };
  is_deeply($deletion->_get_irods_meta_data(),  $expected_meta, 'correct meta data for bam header generated'); 
  
  lives_ok {$deletion->_reload_bam_header($expected_meta)} 'reload bam header back to irods';
  
  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run', 1 unless ($ENV{TOOLS_INSTALLED});
    is_deeply( $loader->_check_meta_data($header), $expected_bam_header_meta, 'correct new bam header meta data' );
  }
  is($loader->rm_file($header), 1, 'delete bam header');
  
  is($deletion->_remove_backups($deletion->bam_file), 1, 'try to remove backup but only one copy');
  is($deletion->_add_backup_removed_meta_data(), 1, 'add backup_removed meta data');
  is($deletion->_remove_irods_file($bam), 1, 'delete copy');
  
  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run', 1 unless ($ENV{TOOLS_INSTALLED});
    throws_ok {$deletion->_generate_bam_header();} qr/Failed to generate bam header: /,
    'cannot generate bam header if bam file does not exist';
  }
}

{
  my $meta = {type      => 'bam',
              id_run    => '1000',
              sample    => [qw(s1 s2)],
              target    => 0,
              'sample_consent_withdrawn_email_sent'=> 1,
             };
  my $loader = new npg_common::irods::Loader(
            file      => 't/data/test.bam',
            collection => $IRODS_TEST_AREA,
            meta_data => $meta,
  );
  $loader->run();

  my $bam = $IRODS_TEST_AREA . '/test.bam';
  my $deletion = npg_common::irods::BamDeletion->new( 
                      bam_file          => $bam,
                      complete_deletion => 1,
                      working_dir       => $working_dir,
                 );

  my $header = $IRODS_TEST_AREA . '/' . $deletion->_header_file_name;

  my $expected_meta = {
          'target' => '0',
          'type' => 'bam_header',
          'sample' => ['s1','s2'],
          'id_run' => ['1000'],
          'sample_consent_withdrawn_email_sent'=> ['1'],
        };
  is_deeply( $deletion->_get_irods_meta_data(), $expected_meta, 'new metadata without rt ticket number');

  my $d;
  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run', 1 unless ($ENV{TOOLS_INSTALLED});
    $d = npg_common::irods::BamDeletion->new( 
                      bam_file          => $bam,
                      complete_deletion => 1,
                      rt_ticket         => 3456,
                      working_dir       => $working_dir,
                   );
    isa_ok($d, 'npg_common::irods::BamDeletion');
    $d->samtools_irods_cmd;
  } 

  SKIP: {
    skip 'samtools_irods executable found', 1 if ($ENV{TOOLS_INSTALLED});
    $d = npg_common::irods::BamDeletion->new( 
                      bam_file          => $bam,
                      complete_deletion => 1,
                      rt_ticket         => 3456,
                      working_dir       => $working_dir,
                      samtools_irods_cmd => $samtools_irods,
                   );
    isa_ok($d, 'npg_common::irods::BamDeletion');
  }

  $expected_meta->{'sample_consent_withdrawn_email_sent'} = 'RT#3456';
  my $bam_header_meta = $d->_get_irods_meta_data();
  is_deeply( $bam_header_meta, $expected_meta, 'new metadata contain rt ticket number');


  lives_ok {$d->_generate_bam_header();} ' bam header generated';
  lives_ok {$d->_reload_bam_header($bam_header_meta);} 'reloaded bam header';

  $expected_bam_header_meta->{'sample_consent_withdrawn_email_sent'} = {'RT#3456' => 1};    

  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run', 1 unless ($ENV{TOOLS_INSTALLED});
    is_deeply( $loader->_check_meta_data( $header ), $expected_bam_header_meta, 'correct new bam header meta data' );
  }

  lives_ok {$d->_remove_irods_file($d->bam_file())} 'removed bam file';

  $d = npg_common::irods::BamDeletion->new( bam_file => $bam,);
  is($d->_remove_backups($d->bam_file()), 1, 'no backup to be removed');
}

{
  my $loader1 = new npg_common::irods::Loader(
     file => 't/data/test.bam', collection => $IRODS_TEST_AREA,);
  $loader1->run();
  my $loader2 = new npg_common::irods::Loader(
     file => 't/data/test.bai', collection => $IRODS_TEST_AREA,);
  $loader2->run();
  
  my $bam = $IRODS_TEST_AREA . '/test.bam';
  my $deletion = npg_common::irods::BamDeletion->new( 
                                                 bam_file => $bam,
                                              );
  ok($deletion->_bai_file_existed(), 'bai file existed');
  
  lives_ok {$deletion->process()} 'remove backups of bam file';
 
  my $deletion2;
  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run', 1 unless ($ENV{TOOLS_INSTALLED});
    $deletion2 = npg_common::irods::BamDeletion->new(
                                     bam_file          => $bam,
                                     complete_deletion => 1,      
                                                    );
    $deletion2->samtools_irods_cmd();
    lives_ok {$deletion2->process()} 'remove bam file completely';
  }
}
   
}; #end of SKIP


sub exist_irods_executables {
   return 0 unless `which ienv`;
   return 0 unless `which imkdir`;
   return 1;
}

sub create_irods_test_area {
  system("imkdir $IRODS_TEST_AREA") == 0 or return 0;
  return 1;
}

END {
  if ($irods_test_area_created) {
    `irm -r $IRODS_TEST_AREA`;    
  }
}

1;
__END__
