use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 4;
use Test::Exception;
use lib qw{t t/lib};
use File::Temp qw(tempdir);
use Cwd;

my $basedir = tempdir( CLEANUP => 1 );
$ENV{TEST_DIR} = $basedir; #so when npg_tracking::illumina::run::folder globs the test directory
$ENV{TEST_FS_RESOURCE} = q{nfs_12}; # as this might be running elsewhere, for now populate from environment variable. In future, write tests which will run if on correct filesystem
use_ok(q{npg_common::role_tests::run_fs_resource});

my $id_run = q{1234};
my $name = q{IL2_1234};
my $run_folder = q{123456_IL2_1234};
my $runfolder_path = qq{$basedir/nfs/sf12/IL2/analysis/123456_IL2_1234};
my $data_subpath = $runfolder_path;

sub delete_staging {
  `rm -rf $basedir/nfs/sf12`;
  return 1;
}

sub create_staging {
  delete_staging();
  `mkdir -p $data_subpath`;
  return 1;
}

my $orig_dir = getcwd();

{
  create_staging();
  my $fs_resource;
  lives_ok  { $fs_resource = npg_common::role_tests::run_fs_resource->new({runfolder_path => $runfolder_path}); } q{created role_test object ok};
  is($fs_resource->runfolder_path(), $runfolder_path, q{runfolder_path found});
  is($fs_resource->fs_resource(), q{nfs_12}, q{correct resource returned});
}

eval { delete_staging(); } or do { carp 'unable to delete staging area'; };
1;
