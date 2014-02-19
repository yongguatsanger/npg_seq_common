#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2009-08-27
# Last Modifi0ed: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use Test::More tests => 42;
use Test::Exception;
use English qw( −no_match_vars );
use Test::Deep;
use File::Temp qw/ tempdir  /;

use_ok('npg_common::irods::Loader');

my $tdir = tempdir( CLEANUP => 1 );
my @path = split '/', $tdir;
my $dir = pop @path;

my $test_base = '/seq/npg/test';
my $test_base1 = '/seq/npg/test1';
my @ids = getpwuid $UID;
my $username = $ids[0];

my $IRODS_TEST_AREA = "$test_base/$dir";
my $IRODS_TEST_AREA1 = "$test_base1/$dir";
my $IRODS_TEST_FILE = "$IRODS_TEST_AREA/RunInfo.xml";
my $TEST_FILE = 't/data/staging/120118_HS22_07385_A_D0GKHACXX/RunInfo.xml';
my $TEST_FILE2 = 't/data/staging/120118_HS22_07385_A_D0GKHACXX/RunInfo2.xml';

my $loader = npg_common::irods::Loader->new(file => 'test/test.file', resource => 'seq');
isa_ok($loader, 'npg_common::irods::Loader', 'object test');

my $EXIST_EXECUTABLES = exist_irods_executables();
my $test_area_created = $EXIST_EXECUTABLES ? create_irods_test_area() : 0;

SKIP: {
  unless ($EXIST_EXECUTABLES) {
    skip 'unable to access iRODS executables', 40;
  }
  elsif (!$test_area_created) { 
    skip 'unable to create iRODS test area (try kinit to log in to iRODS)', 40;
  }

{
  is($loader->new_filename, 'test.file', 'destination file name');
  
  is($loader-> _get_file_md5('t/data/staging/IL19/outgoing/100817_IL19_05169/Data/Intensities/Bustard1.8.1a1_24-08-2010_RTA/GERALD_24-08-2010_RTA/archive/5169_1.bam'), '3ecff85807c7dc73f39f30b3b2a8d6a6', 'correct md5 from corresponding md5 file');
   is($loader-> _get_file_md5('t/data/staging/IL19/outgoing/100817_IL19_05169/Data/Intensities/Bustard1.8.1a1_24-08-2010_RTA/GERALD_24-08-2010_RTA/archive/5169_3.bam'), 'd1e99fdd20aab765f1d8868cb226d02e', 'correct md5 from file itself');
   
  is( $loader->_rm_metadata("test.file", "test_att") , "rmw -d test.file test_att %", "correct rmw imeta command");
  
  is( $loader->_check_add_meta("test.file", "test_att", q{test" value}, {}), q{add -d test.file test_att 'test" value'}, 'correct add command with double quote values');
  is( $loader->_check_add_meta("test.file", "test_att", q{test' value'}, {}), q{add -d test.file test_att "test' value'"}, 'correct add command with single quote values');
  is( $loader->_check_add_meta("test.file", "test_att", q{test" value'}, {}), q{add -d test.file test_att "test' value'"}, 'correct add command with mixed quotation values');
  
    my ($rmw_cmd, $add_cmd) = $loader->_check_add_meta("test.file", "test_att", qq{test" value'}, { test_att => {test => 1}, } );
    is($rmw_cmd, qq{rmw -d test.file test_att %}, 'correct add command with mixed quotation values');
    is($add_cmd, qq{add -d test.file test_att "test' value'"}, 'correct add command with mixed quotation values');

}

{
  my $loader = npg_common::irods::Loader->new(_dry_run=>1, file => 'test/test.file', resource => 'seq');
  $loader->add_meta("test.file", {test_att=>"", test2=>'fred'});
}

{
	my $loader = new npg_common::irods::Loader(file => "$TEST_FILE", new_filename => "$IRODS_TEST_FILE");
	my $r = 0;
	lives_ok{ $r = $loader->run(); } 'no error loading';
	is($r,1, 'return value from the loader is 1 - success');
	
	my $replications_nums = $loader->get_replication_numbers("$IRODS_TEST_FILE");
	my $expected = {rep_nums => [0, 1], obsolete_rep_nums  => [] };
	is_deeply($expected, $replications_nums, 'replications number is correct');
	
	my $obsolete_rep_nums = $loader->_get_rep_num_to_remove("$IRODS_TEST_FILE");
	my $expected_obsolete = [];
	is_deeply($expected_obsolete, $obsolete_rep_nums, 'there are no obsolete replications');

	# we expect it not to load this time
	$loader = new npg_common::irods::Loader(file => "$TEST_FILE", new_filename => "$IRODS_TEST_FILE");
	$r = 0;
	lives_ok{ $r = $loader->run(); } 'no error loading the same file';
	isnt($r,1, 'return value of the loader is not 1 if the same file is loaded');
	
        # try to load a updated version
        $loader = new npg_common::irods::Loader(
                    file         => "$TEST_FILE2",
                    new_filename => "$IRODS_TEST_FILE"
                                              );
	$r = 0;
	lives_ok{ $r = $loader->run(); } 'no error loading an updated file';
	is($r,1, 'return value from the loader is 1 - success');

	$replications_nums = $loader->get_replication_numbers("$IRODS_TEST_FILE");
  # Note that irods is returning [1, 2] or [0, 2] since upgrade, so just check there are two numbers in array
  my $rep_nums = $replications_nums->{rep_nums};
  is(scalar @$rep_nums, 2, q{rep_nums has two entries});
  is($replications_nums->{obsolute_rep_nums}, undef, q{obsolete_rep_nums is empty});

	$obsolete_rep_nums = $loader->_get_rep_num_to_remove("$IRODS_TEST_FILE");
	$expected_obsolete = [];
	is_deeply($expected_obsolete, $obsolete_rep_nums, 'obsolete replications number as expected');
}

{   
	my $loader = new npg_common::irods::Loader(file => "$TEST_FILE", new_filename => "$IRODS_TEST_FILE");
        $loader->remove_one_replication("$IRODS_TEST_FILE", 1);
   
        is($loader->run(), 1, 'one replication successfully removed?');
	is($loader->rm_file("$IRODS_TEST_FILE"),1, 'file removed completely');
}

{
        # Try processing&adding a directory

	my $loader = new npg_common::irods::Loader(file => 't/data/staging/120118_HS22_07385_A_D0GKHACXX/InterOp');
	my $r = 0;
	lives_ok{ $r = $loader->run(); } 'copying interop dir runs ok';
	is($r,1, 'loaded returns 1 - success');

        # then read the resulting collection

	$loader = new npg_common::irods::Loader();
	my $col_list = $loader->get_collection_file_list('InterOp');
	is(scalar keys %$col_list, 7, 'the size of the loaded collection is correct when read back');

        # then delete the collection

	$loader = new npg_common::irods::Loader();
	is($loader->rm_file('InterOp'),1, 'the collection is successfully deleted' );
}

{
        my $file = "$IRODS_TEST_AREA/Changes";
        my $user = $username . '#Sanger1';
        `iput Changes $file`;
        my $other_user = 'public#seq';
        `ichmod read $other_user $file`;

        my $expected =  {'read' => [$other_user], 'own' => [$user]};
        my $loader = npg_common::irods::Loader->new();
        cmp_deeply($loader->get_permissions($file), $expected, 'permissions retrieved');
        my $dry_loader = npg_common::irods::Loader->new(_dry_run => 1);
        lives_ok { $dry_loader->restrict_file_permissions($file)} 'dry run for restricting file permissions lives';
        cmp_deeply($loader->get_permissions($file), $expected, 'permissions not changed');
        lives_ok { $loader->restrict_file_permissions($file)} 'restricting file permissions lives';
        cmp_deeply($loader->get_permissions($file),
          {'own' => [$user]}, 'permissions restricted correctly');

        is($loader->file_exists($file), 1, 'file existence check for existing file');
        my $result;
        lives_ok {$result = $loader->file_exists('/seq/dodo')} 'file existence check does not throw an error';
        ok(!$result, 'file existence check for non-existing file');
}

{
  my $yhuman_bam_filename1 = '10371_1#36.bam';
  my $yhuman_bam_file1 = "t/data/sequence/100818_IL32_10371/Latest_Summary/archive/lane1/$yhuman_bam_filename1";
  my $yhuman_bam_file1_md5 = '3c4881995fa612be0dfe51c4db1b1240';

  my $loader1 = new npg_common::irods::Loader( 
                  file => $yhuman_bam_file1, 
                  collection => $IRODS_TEST_AREA1);

  lives_ok{ my $r = $loader1->run(); } "no error loading $yhuman_bam_file1";
  is($loader1->_get_file_md5($yhuman_bam_file1), $yhuman_bam_file1_md5, 'correct md5 from file itself');
  
  my $other_user = 'public#seq';
  my $file1 = "$IRODS_TEST_AREA1/$yhuman_bam_filename1";
  my $permissions1 = $loader1->get_permissions($file1);
  ok ((grep {$_ eq $other_user} @{$permissions1->{'read'}}), "public read permissions retrieved for $file1");

  my $yhuman_bam_filename2 = '10371_1#36_yhuman.bam';
  my $yhuman_bam_file2 = "t/data/sequence/100818_IL32_10371/Latest_Summary/archive/lane1/$yhuman_bam_filename2";
  my $yhuman_bam_file2_md5 = '385c58cc11ffb6625cd89c65b6d6d213';

  my $loader2 = new npg_common::irods::Loader(
                  file => $yhuman_bam_file2, 
                  collection => $IRODS_TEST_AREA1,
                  chmod_permissions => [q{null public}]);

  lives_ok{ my $r = $loader2->run(); } "no error loading $yhuman_bam_file2 with restricted permissions";
  is($loader2->_get_file_md5($yhuman_bam_file2), $yhuman_bam_file2_md5, 'correct md5 from file itself');
  
  my $file2 = "$IRODS_TEST_AREA1/$yhuman_bam_filename2";
  my $permissions2 = $loader2->get_permissions($file2);
  
  my @plist = $permissions2->{'read'} ? @{$permissions2->{'read'}} : ();
  ok (!(grep {$_ eq $other_user} @{$permissions2->{'read'}}), "no public read permissions retrieved for $file2");
} 
}; #end of skip no irods

sub exist_irods_executables {
  return 0 unless `which ienv`;
  return 0 unless `which imkdir`;
  return 1;
}

sub create_irods_test_area {
  system("imkdir $IRODS_TEST_AREA") == 0 or return 0;
  system("imkdir $IRODS_TEST_AREA1") == 0 or return 0;
  return 1;
}

END {
  if ($test_area_created) {
    eval {system("irm -r $IRODS_TEST_AREA")};
    eval {system("irm -r $IRODS_TEST_AREA1")};
  }
}

1;
__END__
