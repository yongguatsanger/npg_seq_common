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
use Test::More tests => 16;
use Test::Exception;
use Test::Deep;
use File::Temp qw/ tempdir  /;
use Test::MockObject;

use_ok('npg_common::irods::BamConsentWithdrawn');

isa_ok(npg_common::irods::BamConsentWithdrawn->new(), 'npg_common::irods::BamConsentWithdrawn');

is(npg_common::irods::BamConsentWithdrawn->new()->_iquery,
 q{iquest --no-page -z seq "%s/%s" "select COLL_NAME, DATA_NAME where META_DATA_ATTR_NAME = 'sample_consent_withdrawn' and META_DATA_ATTR_VALUE = '1' and DATA_NAME not like '%header.bam%'"}, 'default iquery');
 
is(npg_common::irods::BamConsentWithdrawn->new(zone => undef)->_iquery,
 q{iquest --no-page "%s/%s" "select COLL_NAME, DATA_NAME where META_DATA_ATTR_NAME = 'sample_consent_withdrawn' and META_DATA_ATTR_VALUE = '1' and DATA_NAME not like '%header.bam%'"}, 'iquery without zone');

my $tdir = tempdir( CLEANUP => 1 );
my @path = split '/', $tdir;
my $dname = pop @path;

my $EXIST_EXECUTABLES = exist_irods_executables();

my $IRODS_ENV = $EXIST_EXECUTABLES ? create_irods_test_area(\$dname) : 0;
my $IRODS_TEST_AREA = $dname;
my $IZONE= $IRODS_ENV ? $IRODS_ENV->{'irodsZone'} : undef;
my $username = $IRODS_ENV ? $IRODS_ENV->{'irodsUserName'} : q[]; 
my $imeta = "imeta";
my $ichmod = "ichmod";
my $OTHER_USER = 'public';

SKIP: {

    if ( !$EXIST_EXECUTABLES ) {
        skip 'unable to access iRODS executables', 12;
    } elsif (!$IRODS_ENV) { 
        skip 'unable to create iRODs test area (try kinit to log in to iRODS)', 12;
    }

    my $util = npg_common::irods::BamConsentWithdrawn->new()->_util;
    my @bam_files = qw /1.bam 2.bam 3.bam/;
    my @files = ();
    push @files, @bam_files, qw /1.bai 3.bai 4.bai/;
    foreach my $file (@files) {
      my $source = "$tdir/$file";
      `touch $source`;
      my $target = "$IRODS_TEST_AREA/$file";
      my $iput = "iput";
      `$iput $source $target`;
      $util->file_exists($target) or die "Cannot create $target";
    }

    my $b = npg_common::irods::BamConsentWithdrawn->new(new_bam_files =>
                       ["$IRODS_TEST_AREA/1.bam", "$IRODS_TEST_AREA/2.bam", "$IRODS_TEST_AREA/3.bam"],
                       dry_run => 1,
                                                     );
    is(join(q[ ],@{$b->new_files}),
       "$IRODS_TEST_AREA/1.bam $IRODS_TEST_AREA/1.bai $IRODS_TEST_AREA/2.bam $IRODS_TEST_AREA/3.bam $IRODS_TEST_AREA/3.bai",
       'bai file found');
    lives_ok {$b->process} q[dry run for 'process' when no file has sample_consent_withdrawn_email_sent flag set];

    `$imeta add -d  $IRODS_TEST_AREA/1.bam sample_consent_withdrawn "1"`;
    `$imeta add -d  $IRODS_TEST_AREA/2.bam sample_consent_withdrawn "1"`;
    `$imeta add -d  $IRODS_TEST_AREA/3.bam sample_consent_withdrawn "1"`;
    `$imeta add -d  $IRODS_TEST_AREA/3.bam sample_consent_withdrawn_email_sent "1"`;
     my $other = join q[#], $OTHER_USER, $IZONE;
    `$ichmod write $other $IRODS_TEST_AREA/1.bam`;
    
    my $found = ["$IRODS_TEST_AREA/1.bam", "$IRODS_TEST_AREA/2.bam", "$IRODS_TEST_AREA/3.bam"];

    $b = npg_common::irods::BamConsentWithdrawn->new(dry_run => 1, bam_files => $found);
    is(join(q[ ],@{$b->new_bam_files}), "$IRODS_TEST_AREA/1.bam $IRODS_TEST_AREA/2.bam",
        'bam files with sample_consent_withdrawn_email_sent flag not set found');
    is(join(q[ ],@{$b->new_files}),
        "$IRODS_TEST_AREA/1.bam $IRODS_TEST_AREA/1.bai $IRODS_TEST_AREA/2.bam",
        'full file list');
    lives_ok {$b->_create_rt_ticket} 'dry run for creating an rt ticket';
    lives_ok {$b->process} q[dry run for 'process'];

    my $mock = Test::MockObject->new();
    $mock->fake_new( 'MIME::Lite' );
    $mock->set_true('send');
    my $user = join q[#], $username, $IZONE;

    $b = npg_common::irods::BamConsentWithdrawn->new(zone => undef);
    lives_ok {$b->process} q[live run for 'process'];
    cmp_deeply($b->_util->get_permissions("$IRODS_TEST_AREA/1.bam"),
          {'own' => [$user]}, 'permissions restricted correctly');
    my $get_meta_cmd = "$imeta ls -d $IRODS_TEST_AREA/1.bam";
    my $meta = `$get_meta_cmd`;
    ok($meta =~ /sample_consent_withdrawn_email_sent/, 'sample_consent_withdrawn_email_sent flag set');
    ok(npg_common::irods::BamConsentWithdrawn::_rt_ticket_exists(
                    $util->_check_meta_data("$IRODS_TEST_AREA/1.bam")), 'set flag is recognised');

    $b = npg_common::irods::BamConsentWithdrawn->new(bam_files => $found);
    ok(!@{$b->new_files}, 'no files to process');
    lives_ok {$b->process} q[live run for 'process' where no files found];
};

exit;

sub exist_irods_executables {
   return 0 unless `which ienv`;
   return 0 unless `which imkdir`;
   return 1;
}

sub create_irods_test_area {
  my $dirPtr = shift;

  my $env = `ienv`;
  if (!$env) { return 0; }

  my @lines = split "\n", $env;
  my $envh = {};
  foreach my $line (@lines) {
    if ($line =~ /Version/) {next; }
    my ($key, $value) = $line =~ /\ (\w*)=(.*)$/;
    $envh->{$key} = $value;
  }

  my $dir = ${$dirPtr};
  system("imkdir $dir") == 0 or return 0;
  return $envh;
}

END {
  if ($IRODS_ENV) {
    my @commands = (
      "$imeta rmw -d  $IRODS_TEST_AREA/1.bam sample_consent_withdrawn %",
      "$imeta rmw -d  $IRODS_TEST_AREA/1.bam sample_consent_withdrawn_email_sent %",
      "$imeta rmw -d  $IRODS_TEST_AREA/2.bam sample_consent_withdrawn %",
      "$imeta rmw -d  $IRODS_TEST_AREA/2.bam sample_consent_withdrawn_email_sent %",
      "$imeta rmw -d  $IRODS_TEST_AREA/3.bam sample_consent_withdrawn %",
      "$imeta rmw -d  $IRODS_TEST_AREA/3.bam sample_consent_withdrawn_email_sent %",
      "irm -r $IRODS_TEST_AREA"
                   );
    foreach my $command (@commands) {
      eval {system($command)};
    }
  }
}

1;
