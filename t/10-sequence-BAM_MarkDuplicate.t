#########
# Author:        gq1
# Created:       2009-06-21
#

use strict;
use warnings;
use Test::More tests => 25;
use Test::Exception;

use File::Temp qw(tempdir);
my $temp_dir = tempdir( CLEANUP => 1 );

my $elc_memory_for_deployment = 300;
my $bts_memory_for_deployment = 200;
my $elc_memory_for_production = 3000;
my $bts_memory_for_production = 2000;

use_ok('npg_common::sequence::BAM_MarkDuplicate');

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
            7 unless ($ENV{'TOOLS_INSTALLED'});
  my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam     => 'input.bam',
                 output_bam    => 'output.bam',
                 metrics_json  => 'metrics.json',
                 not_strip_bam_tag => 1,
               );
  isa_ok($bam, 'npg_common::sequence::BAM_MarkDuplicate');
  lives_ok {$bam->temp_dir} 'temp dir generated';
  lives_ok {$bam->metrics_file()} 'temp metrics file';
  lives_ok {$bam->_result} 'result object';
  is($bam->default_java_xmx_elc, $elc_memory_for_production, q{No elc memory supplied so default used});
  is($bam->default_java_xmx_bts, $bts_memory_for_production, q{no bts memory supplied so default used});
  $bam->metrics_file('metrics.txt');
  $bam->temp_dir($temp_dir);
  like($bam->mark_duplicate_cmd(), qr/bammarkduplicates I=input\.bam O=\/dev\/stdout tmpfile=$temp_dir\/ M=metrics\.txt/, 'correct picard command with absolute path to jar');
  is($bam->bamseqchksum_cmd(q{bam}), q{bamseqchksum verbose=0 inputformat=bam > output.bam.seqchksum}, 'correct bamseqchsum command for a bam file');
  is($bam->bamseqchksum_cmd(q{cram}), q{bamseqchksum verbose=0 inputformat=cram > output.cram.seqchksum}, 'correct bamseqchsum command for a cram file with no reference');
       };
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         14 unless ($ENV{'TOOLS_INSTALLED'});
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
               {
                 input_bam     => 't/data/sequence/SecondCall/4392_1.bam',
                 output_bam    => "$temp_dir/output_mk.bam",
                 metrics_json  => "$temp_dir/metrics.json",
                 sort_input    => 1,
                 temp_dir      => $temp_dir,
                 metrics_file  => $temp_dir . '/metrics.txt',
                 not_strip_bam_tag => 1,
                 reference => 't/data/references/E_coli/default/fasta/E-coli-K12.fa',
                 default_java_xmx_elc => $elc_memory_for_deployment,
                 default_java_xmx_bts => $bts_memory_for_deployment,
               });
      my $expected_mark_duplicate_cmd = qq{bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics.txt};
      like($bam->mark_duplicate_cmd(), qr/$expected_mark_duplicate_cmd/, 'correct biobambam command');
      ok( $bam->no_alignment(), 'input bam with alignment');
      $bam->no_alignment(1);
      my $expected_bam_tag_stripper_cmd = qq{INPUT=t/data/sequence/SecondCall/4392_1.bam OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      like($bam->bam_tag_stripper_cmd(), qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');

      $bam->clear_bam_tag_stripper_cmd();  
      $bam->clear_no_alignment();
      $bam->input_bam('t/data/sequence/6062_1#0.bam');
      
      $bam->not_strip_bam_tag(0);
      $bam->no_alignment(0);
       $expected_bam_tag_stripper_cmd = qq{INPUT=/dev/stdin OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      like($bam->bam_tag_stripper_cmd(), qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      # stop qr// in like interpolating READ_NAME_REGEX by enclosing value in \Q..\E
      my $expected_elc_cmd = qq{INPUT=t/data/sequence/6062_1#0.bam OUTPUT=$temp_dir/metrics.txt TMP_DIR=$temp_dir READ_NAME_REGEX='\Q[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\E' VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'};
      like($bam->estimate_library_complexity_cmd(), qr/$expected_elc_cmd/, 'correct elc command');

      my $expected_bamseqchk_cmd = qq{bamseqchksum verbose=0 inputformat=cram reference=t/data/references/E_coli/default/fasta/E-coli-K12.fa > $temp_dir/output_mk.cram.seqchksum};
      is($bam->bamseqchksum_cmd(q{cram}), $expected_bamseqchk_cmd, 'correct bamseqchsum command for a cram file with reference');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      ok ($bam->_result->info->{'Samtools'}, 'samtools version is defined');
      ok ($bam->_result->info->{'Picard-tools'}, 'test picard version is defined');

      lives_ok{$bam->process()} q{Processed OK};
      is (-e "$temp_dir/output_mk.bam", 1, 'BAM file created');      
      is (-e "$temp_dir/output_mk.bai", 1, 'BAM index created');      
      is (-e "$temp_dir/output_mk.bam.md5", 1, 'BAM md5 created');      
      is (-e "$temp_dir/output_mk.flagstat", 1, 'BAM flagstat created');      
      is (-e "$temp_dir/output_mk.cram", 1, 'CRAM file created');
  }    
}

1;
