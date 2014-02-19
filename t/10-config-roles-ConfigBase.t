# $Id$
use strict;
use warnings;
use Carp;
use English qw{-no_match_vars};
use Test::More tests => 14;
use Test::Exception;
use Test::Deep;
use lib qw{t};
use t::util;
use DateTime;

BEGIN {
  use_ok( q{npg_common::config::roles::ConfigBase} );
}

{
  package t::MyConfig;
  use Moose;
  with q{npg_common::config::roles::ConfigBase};

}

local $ENV{HOME} = q{t/data/file_system/home};

{
  my $myconfig;
  lives_ok {
    $myconfig = t::MyConfig->new();
  } q{create object ok};

  isa_ok( $myconfig, q{t::MyConfig}, q{$myconfg} );

  is( $myconfig->conf_name(), q{config.ini}, q{conf_name taken from defaults} );
  is( $myconfig->conf_dir(),  q{data},       q{conf_dir taken from defaults}  );
  is( $myconfig->domain(),    q{live},       q{domain taken from defaults}    );
}

{
  my $myconfig;
  lives_ok {
    $myconfig = t::MyConfig->new( {
      conf_name => q{if_this_is_found_someone_must_be_mad.ini}
    } );
  } q{create object ok};
  

  throws_ok {
    $myconfig->file_config();
  } qr{Searched[ ]for[ ]if_this_is_found_someone_must_be_mad.ini\nIncluded[ ]paths.*\nNo[ ]config[ ]file[ ]found!}, q{cannot find if_this_is_found_someone_must_be_mad.ini};

}

local $ENV{TEST} = q{test};
{
  
  my $myconfig = t::MyConfig->new( {
    conf_name => q{for_ini_test.ini},
    croak_config_read_if_warnings => 1,
  } );

  is( $myconfig->domain(), q{test}, q{domain taken from $ENV{TEST}} );

  my $config_hash;
  lives_ok {
    $config_hash = $myconfig->file_config();
  } q{Obtain config_hash};
  is_deeply( $config_hash, { 'dbhost' => 'ignore_me_please' }, q{section is correct} );

  $myconfig = t::MyConfig->new( {
    conf_name => q{multilevel_config.yml},
    croak_config_read_if_warnings => 1,
    domain => q{live},
  } );

  lives_ok {
    $config_hash = $myconfig->file_config();
  } q{Obtain config_hash};

  is_deeply( $config_hash, {
   'a' => {
     'create_fastq' => 1,
     'create_srf' => 1
   },
   'b' => {
     'create_fastqcheck' => 1,
     'create_gcfreq' => 1
   },
   'c' => {
     'qc_adapter' => 1,
     'qc_gc_fraction' => 1,
     'qc_insert_size' => 1,
     'qc_qX_yield' => 1,
     'qc_sequence_error' => 1,
     'split_nonconsented_sequence' => 1
   },
   'd' => {
     'archive_to_irods' => 1,
     'upload_auto_qc_to_qc_database' => 1,
     'upload_fastqcheck_to_qc_database' => 1,
     'upload_illumina_analysis_to_qc_database' => 1,
     'upload_summary_to_qc_database' => 1
   },
   'e' => {
     'harold_calibration_tables' => 1,
     'make_gerald_tiles' => 1
   },
   'f' => {
     'harold_recalibration' => 1,
     'make_gerald_qtable' => 1
   },
   'g' => {
     'bam_generation' => 1,
     'qc_contamination' => 1
   },
   'h' => {
     'bam_markduplicate' => 1,
     'create_md5' => 1
   }
   }, q{section is correct} );

  $myconfig = t::MyConfig->new( {
    conf_name => q{for_test_no_test_section.ini},
    croak_config_read_if_warnings => 1,
  } );

  throws_ok {
    $config_hash = $myconfig->file_config();
  } qr{Config[ ]warnings[ ]found[ ]-[ ]croaking:\n\ttest[ ]section[ ]not[ ]found[.][ ]Looking[ ]in[ ]variants[.]\n}, q{no test section and croak on warnings};


}

1;