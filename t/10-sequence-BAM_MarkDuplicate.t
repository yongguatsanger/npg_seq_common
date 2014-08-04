#########
# Author:        gq1
# Created:       2009-06-21
#

use strict;
use warnings;
use Test::More tests => 114;
use Test::Exception;
use Test::Deep;

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
            9 unless ($ENV{'TOOLS_INSTALLED'});
  my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam     => 'input.bam',
                 output_bam    => 'output.bam',
                 metrics_json  => 'metrics.json',
                 not_strip_bam_tag => 1,
                 no_alignment => 0,
                 human_split => 'all', 
               );
  isa_ok($bam, 'npg_common::sequence::BAM_MarkDuplicate');
  lives_ok {$bam->temp_dir} 'temp dir generated';
  lives_ok {$bam->metrics_file()} 'temp metrics file';
  lives_ok {$bam->_result} 'result object';
  is($bam->default_java_xmx_elc, $elc_memory_for_production, q{no elc memory supplied so default used});
  is($bam->default_java_xmx_bts, $bts_memory_for_production, q{no bts memory supplied so default used});
  $bam->metrics_file('metrics.txt');
  $bam->temp_dir($temp_dir);
  like($bam->mark_duplicate_cmd(), qr/bammarkduplicates I=input\.bam O=\/dev\/stdout tmpfile=$temp_dir\/ M=metrics\.txt/, 'correct picard command with absolute path to jar');
  is($bam->bamseqchksum_cmd(q{bam}), q{bamseqchksum verbose=1 inputformat=bam}, 'correct bamseqchksum command for a bam file');
  is($bam->bamseqchksum_cmd(q{cram}), q{bamseqchksum verbose=1 inputformat=cram}, 'correct bamseqchksum command for a cram file with no reference');
  
       };
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         36 unless ($ENV{'TOOLS_INSTALLED'});
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
               {
                 input_bam     => 't/data/sequence/SecondCall/4392_1.bam',
                 output_bam    => "$temp_dir/output_mk.bam",
                 metrics_json  => "$temp_dir/metrics.json",
                 sort_input    => 1,
                 temp_dir      => $temp_dir,
                 metrics_file  => $temp_dir . '/metrics.txt',
                 not_strip_bam_tag => 1,
                 reference => 't/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta',
                 default_java_xmx_elc => $elc_memory_for_deployment,
                 default_java_xmx_bts => $bts_memory_for_deployment,
                 human_split => 'all', 
               });
      my $expected_mark_duplicate_cmd = qq{bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics.txt};
      like($bam->mark_duplicate_cmd(), qr/$expected_mark_duplicate_cmd/, 'correct biobambam command');
      ok($bam->no_alignment(), 'input bam with alignment');
      $bam->no_alignment(1);
      my $expected_bam_tag_stripper_cmd = qq{INPUT=t/data/sequence/SecondCall/4392_1.bam OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      like($bam->bam_tag_stripper_cmd(), qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');

      $bam->clear_bam_tag_stripper_cmd();  
      $bam->clear_no_alignment();
      $bam->input_bam('t/data/sequence/plasmodium.bam');
      
      $bam->not_strip_bam_tag(0);
      $bam->no_alignment(0);
       $expected_bam_tag_stripper_cmd = qq{INPUT=/dev/stdin OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      my $bam_tag_stripper_cmd = $bam->bam_tag_stripper_cmd();
      like($bam_tag_stripper_cmd, qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      # stop qr// in like interpolating READ_NAME_REGEX by enclosing value in \Q..\E
      my $expected_elc_cmd = qq{INPUT=t/data/sequence/plasmodium.bam OUTPUT=$temp_dir/metrics.txt TMP_DIR=$temp_dir READ_NAME_REGEX='\Q[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\E' VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'};
      like($bam->estimate_library_complexity_cmd(), qr/$expected_elc_cmd/, 'correct elc command');

      my $bam_pb_cal_cmd = $bam->pb_cal_cmd();
      is($bam_pb_cal_cmd, q{/software/solexa/pkg/pb_calibration/v10.15/bin/calibration_pu}, 'correct pb_cal comand');

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file, using newer samtools than current');

      my $expected_bamseqchk_cmd = qq{bamseqchksum verbose=1 inputformat=cram reference=t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta};
      is($bam->bamseqchksum_cmd(q{cram}), $expected_bamseqchk_cmd, 'correct bamseqchksum command for a cram file with reference');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      ok ($bam->_result->info->{'Samtools'}, 'samtools version is defined for an unaligned bam');
      ok ($bam->_result->info->{'Picard-tools'}, 'test picard version is defined for an unaligned bam');

      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;

      my $expected_tee_cmd = qq{set -o pipefail;/software/hpag/biobambam/0.0.147/bin/bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics.txt level='0' | $bam_tag_stripper_cmd | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.md5.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.flagstat.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.bamcheck.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.bschk.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.index.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.pb_cal.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.scramble.fifo};
      $expected_tee_cmd .= qq{ > $temp_dir/output_mk.bam};

      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly');

      my $cram_seqchksum_file_name_mk = qq{$temp_dir/output_mk.cram.seqchksum};
      my $cram_seqchksum_fifo_name_mk = qq{$temp_dir/output_mk.cram.seqchksum.fifo};
      my $cram_file_name_mk = qq{$temp_dir/output_mk.cram};
      my $cram_fifo_name_mk = qq{$temp_dir/output_mk.cram.fifo};
      my $bam_file_name_mk = qq{$temp_dir/output_mk.bam};
      my $bam_seqchksum_file_name_mk = qq{$temp_dir/output_mk.bam.seqchksum};
      my $bam_seqchksum_fifo_name_mk = qq{$temp_dir/output_mk.bam.seqchksum.fifo};

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-"}. qq{ > $temp_dir/output_mk.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_mk.flagstat};

      my $expected_bamcheck_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.bamcheck.fifo | };
      $expected_bamcheck_cmd .= qq{$bam_bamcheck_cmd > $temp_dir/output_mk.bamcheck};

      my $expected_index_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.index.fifo | /software/solexa/pkg/samtools/samtools-$samtools_version/samtools index /dev/stdin /dev/stdout > $temp_dir/output_mk.bai};

      my $expected_pb_cal_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.pb_cal.fifo | };
      $expected_pb_cal_cmd .= qq{$bam_pb_cal_cmd -p t/data/sequence/plasmodium -filter-bad-tiles 2 -};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{bamseqchksum verbose=1 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ | tee $temp_dir/output_mk.bam.seqchksum.fifo > $temp_dir/output_mk.bam.seqchksum};
  
      my $expected_scramble_cmd = qq{/software/badger/bin/scramble -I bam -O cram < $temp_dir/output_mk.bam.scramble.fifo -r t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta };
      $expected_scramble_cmd .= qq{| tee $temp_dir/output_mk.cram.fifo > $temp_dir/output_mk.cram};

      my $expected_cramchksum_cmd =  qq{set -o pipefail; cat $cram_fifo_name_mk | bamseqchksum verbose=1 inputformat=cram reference=t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta };
      $expected_cramchksum_cmd .= qq{| tee $cram_seqchksum_fifo_name_mk > $cram_seqchksum_file_name_mk};

      my $expected_diff_cmd = qq{diff $bam_seqchksum_fifo_name_mk $cram_seqchksum_fifo_name_mk};

      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_scramble_cmd;
      push  @expected_fork_cmds, $expected_cramchksum_cmd;
      push  @expected_fork_cmds, $expected_diff_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_bamcheck_cmd;
      push  @expected_fork_cmds, $expected_index_cmd;
      push  @expected_fork_cmds, $expected_pb_cal_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for fork generated correctly');

      lives_ok{$bam->process()} q{Processed OK};

      is (-e "$temp_dir/output_mk.bam.md5.fifo", 1, 'md5 FIFO created');
      is (-e "$temp_dir/output_mk.bam.flagstat.fifo", 1, 'flagstat FIFO created');
      is (-e "$temp_dir/output_mk.bam.bamcheck.fifo", 1, 'bamcheck FIFO created');
      is (-e "$temp_dir/output_mk.bam.index.fifo", 1, 'index FIFO created');
      is (-e "$temp_dir/output_mk.bam.pb_cal.fifo", 1, 'pb_cal FIFO created');
      is (-e "$temp_dir/output_mk.bam.scramble.fifo", 1, 'scramble FIFO created');
      is (-e "$temp_dir/output_mk.cram.fifo", 1, 'cram FIFO created');
      is (-e "$temp_dir/output_mk.bam.bschk.fifo", 1, 'bamseqchksum input FIFO created');
      is (-e "$temp_dir/output_mk.bam.seqchksum.fifo", 1, 'bamseqchksum output FIFO created');
      is (-e "$temp_dir/output_mk.cram.seqchksum.fifo", 1, 'bamseqchksum output FIFO created');

      is (!-z "$temp_dir/output_mk.bam", 1, 'BAM file created with contents');      
      is (!-z "$temp_dir/output_mk.bai", 1, 'BAM index created with contents');      
      is (!-z "$temp_dir/output_mk.bam.md5", 1, 'BAM md5 created with contents');      
      is (!-z "$temp_dir/output_mk.flagstat", 1, 'BAM flagstat created with contents');      
      is (!-z "$temp_dir/output_mk_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents');
      is (!-z "$temp_dir/output_mk_quality_cycle_surv.txt", 1, 'Quality surv created with contents');
      is (!-z "$temp_dir/output_mk_quality_error.txt", 1, 'Quality error table created with contents');
      is (!-z "$temp_dir/output_mk.cram", 1, 'CRAM file created with contents');
      is (-e "$temp_dir/output_mk.cram.fifo", 1, 'CRAM fifo created');
      is (-e "$temp_dir/output_mk.cram.seqchksum.fifo", 1, 'CRAM seqchksum fifo created');
      is (!-z "$temp_dir/output_mk.bam.seqchksum", 1, 'BAM seqchksum file created with contents');
      is (!-z "$temp_dir/output_mk.cram.seqchksum", 1, 'CRAM seqchksum file created with contents');
  }    
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         32 unless ($ENV{'TOOLS_INSTALLED'});
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
               {
                 input_bam     => 't/data/sequence/unaligned.bam',
                 output_bam    => "$temp_dir/output_no_align.bam",
                 metrics_json  => "$temp_dir/metrics_no_align.json",
                 sort_input    => 1,
                 temp_dir      => $temp_dir,
                 metrics_file  => $temp_dir . '/metrics_no_align.txt',
                 not_strip_bam_tag => 0,
                 default_java_xmx_elc => $elc_memory_for_deployment,
                 default_java_xmx_bts => $bts_memory_for_deployment,
                 no_alignment => 1,
               });
      my $expected_mark_duplicate_cmd = qq{bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics_no_align.txt};
      like($bam->mark_duplicate_cmd(), qr/$expected_mark_duplicate_cmd/, 'correct biobambam command');
      ok( $bam->no_alignment(), 'input bam without alignment');
      my $expected_bam_tag_stripper_cmd = qq{INPUT=t/data/sequence/unaligned.bam OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      my $bam_tag_stripper_cmd = $bam->bam_tag_stripper_cmd();
      like($bam_tag_stripper_cmd, qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      # stop qr// in like interpolating READ_NAME_REGEX by enclosing value in \Q..\E
      my $expected_elc_cmd = qq{INPUT=t/data/sequence/unaligned.bam OUTPUT=$temp_dir/metrics_no_align.txt TMP_DIR=$temp_dir READ_NAME_REGEX='\Q[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\E' VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'};
      like($bam->estimate_library_complexity_cmd(), qr/$expected_elc_cmd/, 'correct elc command');

      is($bam->bamseqchksum_cmd(q{bam}), qq{bamseqchksum verbose=1 inputformat=bam}, 'correct bamseqchksum command for a bam file with reference but no alignment');
      is($bam->bamseqchksum_cmd(q{cram}), qq{bamseqchksum verbose=1 inputformat=cram}, 'correct bamseqchksum command for a cram file with reference but no alignment');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;

      is ($bam->_result->info->{'Picard-tools'}, undef, 'test picard version is not defined if no_alignment flag used');

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file');

      my $expected_tee_cmd = qq{set -o pipefail;$bam_tag_stripper_cmd | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_no_align.bam.md5.fifo $temp_dir/output_no_align.bam.flagstat.fifo $temp_dir/output_no_align.bam.bamcheck.fifo $temp_dir/output_no_align.bam.bschk.fifo > $temp_dir/output_no_align.bam};
      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly if no_alignment flag used');

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-"}. qq{ > $temp_dir/output_no_align.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_no_align.flagstat};

      my $expected_bamcheck_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.bamcheck.fifo | };
      $expected_bamcheck_cmd .= qq{$bam_bamcheck_cmd > $temp_dir/output_no_align.bamcheck};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{bamseqchksum verbose=1 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ > $temp_dir/output_no_align.bam.seqchksum};
  
      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_bamcheck_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for ForkManager generated correctly');

      lives_ok{$bam->process()} q{Processed OK};
      is (-e "$temp_dir/output_no_align.bam.md5.fifo", 1, 'md5 FIFO created');
      is (-e "$temp_dir/output_no_align.bam.flagstat.fifo", 1, 'flagstat FIFO created');
      is (-e "$temp_dir/output_no_align.bam.bamcheck.fifo", 1, 'bamcheck FIFO created');
      is (!-e "$temp_dir/output_no_align.bam.index.fifo", 1, 'index FIFO NOT created');
      is (!-e "$temp_dir/output_no_align.bam.pb_cal.fifo", 1, 'pb_cal FIFO NOT created');
      is (!-e "$temp_dir/output_no_align.bam.scramble.fifo", 1, 'scramble FIFO NOT created');
      is (-e "$temp_dir/output_no_align.bam.bschk.fifo", 1, 'bamseqchksum input FIFO created');
      is (!-e "$temp_dir/output_no_align.bam.seqchksum.fifo", 1, 'bamseqchksum output FIFO NOT created');
      is (!-e "$temp_dir/output_no_align.cram.seqchksum.fifo", 1, 'bamseqchksum output FIFO NOT created');

      is (!-z "$temp_dir/output_no_align.bam", 1, 'BAM file created with contents');      
      is (!-e "$temp_dir/output_no_align.bai", 1, 'BAM index NOT created if no_alignment flag used');      
      is (!-z "$temp_dir/output_no_align.bam.md5", 1, 'BAM md5 created with contents');      
      is (!-z "$temp_dir/output_no_align.flagstat", 1, 'BAM flagstat created with contents');      
      is (!-z "$temp_dir/output_no_align_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents');
      is (!-z "$temp_dir/output_no_align_quality_cycle_surv.txt", 1, 'Quality surv created with contents');
      is (!-z "$temp_dir/output_no_align_quality_error.txt", 1, 'Quality error table created with contents');
      is (!-e "$temp_dir/output_no_align.cram", 1, 'CRAM file NOT created if no_alignment flag used');
      is (!-e "$temp_dir/output_no_align.cram.seqchksum", 1, 'CRAM seqchksum file NOT created if no_alignment flag used');
      is (!-z "$temp_dir/output_no_align.bam.seqchksum", 1, 'BAM seqchksum file created with contents if no_alignment flag used');
  }
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         33 unless ($ENV{'TOOLS_INSTALLED'});
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
               {
                 input_bam     => 't/data/sequence/phix.bam',
                 output_bam    => "$temp_dir/output_phix.bam",
                 metrics_json  => "$temp_dir/metrics_phix.json",
                 sort_input    => 1,
                 temp_dir      => $temp_dir,
                 metrics_file  => $temp_dir . '/metrics_phix.txt',
                 not_strip_bam_tag => 0,
                 default_java_xmx_elc => $elc_memory_for_deployment,
                 default_java_xmx_bts => $bts_memory_for_deployment,
                 human_split => 'phix', 
               });
      my $expected_mark_duplicate_cmd = qq{bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics_phix.txt};
      like($bam->mark_duplicate_cmd(), qr/$expected_mark_duplicate_cmd/, 'correct biobambam command');
      ok(!$bam->no_alignment(), 'input PhiX bam with alignment');
      
      my $expected_bam_tag_stripper_cmd = qq{INPUT=/dev/stdin OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      my $bam_tag_stripper_cmd = $bam->bam_tag_stripper_cmd();
      like($bam_tag_stripper_cmd, qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      is($bam->bamseqchksum_cmd(q{bam}), qq{bamseqchksum verbose=1 inputformat=bam}, 'correct bamseqchksum command for a bam file with reference but no alignment');
      is($bam->bamseqchksum_cmd(q{cram}), qq{bamseqchksum verbose=1 inputformat=cram}, 'correct bamseqchksum command for a cram file with reference but no alignment');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file');

      my $expected_tee_cmd = qq{set -o pipefail;/software/hpag/biobambam/0.0.147/bin/bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics_phix.txt level='0' | $bam_tag_stripper_cmd | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_phix.bam.md5.fifo $temp_dir/output_phix.bam.flagstat.fifo $temp_dir/output_phix.bam.bamcheck.fifo $temp_dir/output_phix.bam.bschk.fifo $temp_dir/output_phix.bam.index.fifo $temp_dir/output_phix.bam.pb_cal.fifo > $temp_dir/output_phix.bam};
      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly for PhiX');

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-"}. qq{ > $temp_dir/output_phix.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_phix.flagstat};

      my $expected_bamcheck_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.bamcheck.fifo | };
      $expected_bamcheck_cmd .= qq{$bam_bamcheck_cmd > $temp_dir/output_phix.bamcheck};

      my $expected_index_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.index.fifo | /software/solexa/pkg/samtools/samtools-$samtools_version/samtools index /dev/stdin /dev/stdout > $temp_dir/output_phix.bai};

      my $bam_pb_cal_cmd = $bam->pb_cal_cmd();
      my $expected_pb_cal_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.pb_cal.fifo | };
      $expected_pb_cal_cmd .= qq{$bam_pb_cal_cmd -p t/data/sequence/phix -filter-bad-tiles 2 -};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{bamseqchksum verbose=1 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ > $temp_dir/output_phix.bam.seqchksum};
 
      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_bamcheck_cmd;
      push  @expected_fork_cmds, $expected_index_cmd;
      push  @expected_fork_cmds, $expected_pb_cal_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for ForkManager generated correctly');

      lives_ok{$bam->process()} q{Processed OK};

      is (-e "$temp_dir/output_phix.bam.md5.fifo", 1, 'md5 FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.flagstat.fifo", 1, 'flagstat FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.bamcheck.fifo", 1, 'bamcheck FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.index.fifo", 1, 'index FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.pb_cal.fifo", 1, 'pb_cal FIFO created for PhiX');
      is (!-e "$temp_dir/output_phix.bam.scramble.fifo", 1, 'scramble FIFO NOT created for PhiX');
      is (-e "$temp_dir/output_phix.bam.bschk.fifo", 1, 'bamseqchksum input FIFO created for PhiX');
      is (!-e "$temp_dir/output_phix.bam.seqchksum.fifo", 1, 'bamseqchksum output FIFO NOT created for PhiX');
      is (!-e "$temp_dir/output_phix.cram.seqchksum.fifo", 1, 'bamseqchksum output FIFO NOT created for PhiX');

      is (!-z "$temp_dir/output_phix.bam", 1, 'BAM file created with contents for PhiX');      
      is (!-z "$temp_dir/output_phix.bai", 1, 'BAM index created with contents for PhiX');      
      is (!-z "$temp_dir/metrics_phix.bam.json", 1, 'metrics json created with contents for PhiX');      
      is (!-z "$temp_dir/metrics_phix.txt", 1, 'metrics txt created with contents for PhiX');      
      is (!-z "$temp_dir/output_phix.bam.md5", 1, 'BAM md5 created with contents for PhiX');      
      is (!-z "$temp_dir/output_phix.flagstat", 1, 'BAM flagstat created with contents for PhiX');      
      is (!-z "$temp_dir/13388_2#40_phix_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents for PhiX');
      is (!-z "$temp_dir/13388_2#40_phix_quality_cycle_surv.txt", 1, 'Quality surv created with contents for PhiX');
      is (!-z "$temp_dir/13388_2#40_phix_quality_error.txt", 1, 'Quality error table created with contents for PhiX');

      is (!-e "$temp_dir/output_no_align.cram", 1, 'CRAM file NOT created for PhiX');
      is (!-e "$temp_dir/output_no_align.cram.seqchksum", 1, 'CRAM seqchksum file NOT created for PhiX');
      is (!-z "$temp_dir/output_no_align.bam.seqchksum", 1, 'BAM seqchksum file created with contents for PhiX');
  }
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         3 unless ($ENV{'TOOLS_INSTALLED'});
    my $incorrect_cmd = 'not_md5sum -b';
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
               {
                 input_bam     => 't/data/sequence/plasmodium.bam',
                 output_bam    => "$temp_dir/output_plasmodium.bam",
                 metrics_json  => "$temp_dir/metrics_plasmodium.json",
                 sort_input    => 1,
                 temp_dir      => $temp_dir,
                 metrics_file  => $temp_dir . '/metrics_plasmodium.txt',
                 not_strip_bam_tag => 0,
                 default_java_xmx_elc => $elc_memory_for_deployment,
                 default_java_xmx_bts => $bts_memory_for_deployment,
                 create_md5_cmd => $incorrect_cmd, 
               });
            lives_ok{$bam} q{Object created with incorrect md5 command};
            is($bam->create_md5_cmd(), $incorrect_cmd, 'Incorrect command created OK');
            throws_ok{$bam->process()} qr/exit/, q{Processing exits when one child process exits because a command inside the pipe is incorrect};
     }
}
 
1;
