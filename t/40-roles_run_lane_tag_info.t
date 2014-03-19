use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 22;
use Test::Exception;
use Test::Deep;
use lib qw{t/lib};
use Cwd;
use File::Temp qw{ tempdir };

BEGIN {
  use_ok( q{npg_common::roles::run::lane::tag_info} );
  use_ok( q{npg_common::role_tests::run_lane_tag_info} );
}

my $basedir = tempdir(CLEANUP => 1);
$ENV{TEST_DIR} = $basedir; #so when npg_tracking::illumina::run::folder globs the test directory
my $id_run = q{1234};
my $name = q{IL2_1234};
my $run_folder = q{123456_IL2_1234};
my $runfolder_path = qq{$basedir/staging/IL2/analysis/123456_IL2_1234};
my $config_path = $runfolder_path . q{/Config};

sub delete_staging {
  `rm -rf $basedir/staging`;
  return 1;
}

sub create_staging {
  delete_staging();
  `mkdir -p $config_path`;
  `cp t/data/recipes/Recipe_GA2-PEM_MP_2x76Cycle+8_v7.7.xml $runfolder_path/`;
  `cp t/data/recipes/TileLayout.xml $config_path/`;
  `cp t/data/recipes/lane_tag_files/lane_1.tag $runfolder_path/`;
  return 1;
}


my $orig_dir = getcwd();

{
  my $lti;
  lives_ok {
    $lti = npg_common::role_tests::run_lane_tag_info->new();
  } q{create object ok};
  isa_ok( $lti, q{npg_common::role_tests::run_lane_tag_info}, q{$lti} );

  throws_ok {
    $lti->add_lane_tag_pair( q{String} );
  } qr{(does|did)[ ]not[ ]pass[ ](its|container)[ ]type[ ]constraint[ ](.+failed.+)?'HashRef[[]Int[]]'}, q{can't push non_hashref via $lti->add_lane_tag_pair()};

  is( $lti->count_lane_tag_pairs(), 0, q{$lti->count_lane_tag_pairs() is 0} );

  throws_ok {
    $lti->add_lane_tag_pair( { position => q{String} } );
  } qr{(does|did)[ ]not[ ]pass[ ](its|container)[ ]type[ ]constraint[ ](.+failed.+)?'HashRef[[]Int[]]'}, q{can't push hashref with non-integer values via $lti->add_lane_tag_pair()};

  ok( $lti->has_no_lane_tag_pairs(), q{$lti->has_no_lane_tag_pairs() is true} );

  throws_ok {
    $lti->add_lane_tag_pairs();
  } qr{No[ ]lane[ ]provided[.][ ]Please[ ]supply[ ]a[ ]lane[ ]position,[ ]or[ ]use[ ]add_all_lane_tag_pairs}, q{no lane position provided};

  lives_ok {
    $lti->lane( 1 );
  } q{no croak with a lane position provided via the lane accessor};

}

{
  my $lti;

  create_staging();

  lives_ok {
    $lti = npg_common::role_tests::run_lane_tag_info->new({
      run_folder => $run_folder,
      runfolder_path => $runfolder_path,
    });
  } q{create object ok};

  ok( $lti->is_multiplexed_lane( 1 ), q{lane 1 is multiplexed} );
  ok( ! $lti->is_multiplexed_lane( 4 ), q{lane 4 is not multiplexed} );

  lives_ok {
    $lti->add_all_lane_tag_pairs();
  } q{$lti->add_all_lane_tag_pairs() ok};

  is( $lti->count_lane_tag_pairs(), 7, q{7 lane_tag pairs found} );

  my @pair_array = $lti->all_lane_tag_pairs();
  my $test_pair_array = [
    { position => 1, tag_index => 0},
    { position => 1, tag_index => 1},
    { position => 1, tag_index => 2},
    { position => 1, tag_index => 3},
    { position => 1, tag_index => 4},
    { position => 1, tag_index => 5},
    { position => 1, tag_index => 6},
  ];

  is_deeply( \@pair_array, $test_pair_array, q{pair data array as expected} );


  is(
    $lti->create_array_string( 1000,1001,1002,1003,1006,3000,3300,3301,3302,3303,3304,3305,3306,3998,3999,4000,4001,4002),
    q{[1000-1003,1006,3000,3300-3306,3998-4002]},
    q{_create_array_string does what it should},
  );

  is( $lti->lsf_job_array_from_lane_tag_pairs(), q{[1000-1006]}, q{correct lsf_job_array obtained for this run} );

  is( $lti->position_decode_string(),  q{`echo $LSB_JOBINDEX/1000 | bc`}, q{position_decode_string is correct}  );
  is( $lti->tag_index_decode_string(), q{`echo $LSB_JOBINDEX%1000 | bc`}, q{tag_index_decode_string is correct} );
}
{

  my $lti;

  create_staging();
  `cp t/data/recipes/lane_tag_files/lane_8_with_too_big_index_number.tag $runfolder_path/lane_8.tag`;

  lives_ok {
    $lti = npg_common::role_tests::run_lane_tag_info->new({
      run_folder => $run_folder,
      runfolder_path => $runfolder_path,
    });
  } q{create object ok};

  throws_ok { $lti->add_all_lane_tag_pairs() } qr/failed[ ](for[ ].*)?with value 6000/, 
              q{croak as an index number is too big};
}

eval { delete_staging(); } or do { carp 'unable to delete staging area'; };

1;
