use strict;
use warnings;
use Test::More tests => 9;
use Test::Exception;
use Perl6::Slurp;
use lib qw{t/lib};

BEGIN {
  use_ok(q{npg_common::role_tests::log});
}
my $log_file_name = q{/tmp/test_log_file};
{
  my $logger;
  lives_ok { $logger = npg_common::role_tests::log->new({}); } q{object created ok};
  lives_ok { $logger->log(q{test}); } q{no problem just doing a log with no memory logging and no log to a file};
  is($logger->memory_log(), q{}, q{as require_memory_log not set, no log in memory});
}
{
  my $logger = npg_common::role_tests::log->new({
    require_memory_log => 1,
  });
  lives_ok { $logger->log(q{test}); } q{no problem just doing a log with memory logging and no log to a file};
  like($logger->memory_log(), qr{test\n\z}, q{memory_log set ok});
}
{
  my $logger = npg_common::role_tests::log->new({
    log_file_name => q{test_log_file},
    log_file_path => q{/tmp},
  });
  lives_ok { $logger->log(q{test}); } q{no problem logging to a file};
  ok((-e $log_file_name), q{log file exists});
  my $log_file = slurp $log_file_name;
  like($log_file, qr{\A\[\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\][ ]test\n\z},
   q{log file contains expected log});
}

END {
  if (-e $log_file_name) {
    qx{rm -f /tmp/test_log_file};
  } 
}
1;