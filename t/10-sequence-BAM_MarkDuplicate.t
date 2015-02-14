#########
# Author:        gq1
# Created:       2009-06-21
#
use strict;
use warnings;
use Test::More tests => 122;
use Test::Exception;
use Test::Deep;
use Cwd;

use File::Temp qw(tempdir);
my $temp_dir = tempdir( CLEANUP => 1 );


my $elc_memory_for_deployment = 300;
my $bts_memory_for_deployment = 200;
my $elc_memory_for_production = 16000;
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
  like($bam->mark_duplicate_cmd(), qr/bammarkduplicates2 I=input\.bam O=\/dev\/stdout tmpfile=$temp_dir\/ M=metrics\.txt/, 'correct picard command with absolute path to jar');
  like($bam->bamseqchksum_cmd(q{bam}), qr{\Qbamseqchksum verbose=0 inputformat=bam\E}, 'correct bamseqchksum command for a bam file');
  like($bam->bamseqchksum_cmd(q{cram}), qr{\Qbamseqchksum verbose=0 inputformat=cram\E}, 'correct bamseqchksum command for a cram file with no reference');
  
       };
}

my $bammarkduplicates = `which bammarkduplicates2`;
$bammarkduplicates = `readlink -f $bammarkduplicates`;
chomp $bammarkduplicates;
my $bamseqchksum = `which bamseqchksum`;
$bamseqchksum = `readlink -f $bamseqchksum`;
chomp $bamseqchksum;
my $scramble = `which scramble`;
$scramble = `readlink -f $scramble`;
chomp $scramble;
my $cram_index = `which cram_index`;
$cram_index = `readlink -f $cram_index`;
chomp $cram_index;
{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         40 unless ($ENV{'TOOLS_INSTALLED'});
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
                 replace_file => 1,
               });
      my $expected_mark_duplicate_cmd = qq{bammarkduplicates2 I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics.txt};
      like($bam->mark_duplicate_cmd(), qr/$expected_mark_duplicate_cmd/, 'correct biobambam command');
      ok($bam->no_alignment(), 'input bam with alignment');
      $bam->no_alignment(1);
      my $expected_bam_tag_stripper_cmd = qq{INPUT=t/data/sequence/SecondCall/4392_1.bam OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      like($bam->bam_tag_stripper_cmd(), qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');

      $bam->clear_bam_tag_stripper_cmd();  
      $bam->clear_no_alignment();
      $bam->not_strip_bam_tag(0);
      $bam->no_alignment(0); ## has alignments
      
      $bam->input_bam('t/data/sequence/5551_3#6.bam'); 
      $expected_bam_tag_stripper_cmd = qq{INPUT=/dev/stdin OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      like($bam->bam_tag_stripper_cmd(), qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      # stop qr// in like interpolating READ_NAME_REGEX by enclosing value in \Q..\E
      my $expected_elc_cmd = qq{INPUT=t/data/sequence/5551_3#6.bam OUTPUT=$temp_dir/metrics.txt TMP_DIR=$temp_dir READ_NAME_REGEX='\Q[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\E' VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'};
      like($bam->estimate_library_complexity_cmd(), qr/$expected_elc_cmd/, 'correct elc command');

      $bam->clear_bam_tag_stripper_cmd();  
      my $current_dir = getcwd();
      system "cp -pv $current_dir/t/data/sequence/plasmodium.bam $temp_dir";
      $bam->input_bam("$temp_dir/plasmodium.bam");
      
      $expected_bam_tag_stripper_cmd = qq{INPUT=/dev/stdin OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      my $bam_tag_stripper_cmd = $bam->bam_tag_stripper_cmd();
      like($bam_tag_stripper_cmd, qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      # stop qr// in like interpolating READ_NAME_REGEX by enclosing value in \Q..\E
      $expected_elc_cmd = qq{INPUT=$temp_dir/plasmodium.bam OUTPUT=$temp_dir/metrics.txt TMP_DIR=$temp_dir READ_NAME_REGEX='\Q[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\E' VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'};
      like($bam->estimate_library_complexity_cmd(), qr/$expected_elc_cmd/, 'correct elc command');

      my $bam_pb_cal_cmd = $bam->pb_cal_cmd();
      like($bam_pb_cal_cmd, qr{/software/solexa/pkg/pb_calibration/\S+/bin/calibration_pu}, 'correct pb_cal comand');

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file, using newer samtools than current');

      my $expected_bamseqchk_cmd = qq{$bamseqchksum verbose=0 inputformat=cram reference=t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta};
      is($bam->bamseqchksum_cmd(q{cram}), $expected_bamseqchk_cmd, 'correct bamseqchksum command for a cram file with reference');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      ok ($bam->_result->info->{'Samtools'}, 'samtools version is defined for an unaligned bam');
      ok ($bam->_result->info->{'Picard-tools'}, 'test picard version is defined for an unaligned bam');

      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;

      my $expected_tee_cmd = qq{set -o pipefail;$bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics.txt level=0 | $bam_tag_stripper_cmd | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.md5.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.flagstat.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.bamcheck.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.bschk.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.alt.bschk.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.index.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.pb_cal.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.scramble.fifo};
      $expected_tee_cmd .= qq{ > $temp_dir/output_mk.bam};

      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly');

      my $cram_seqchksum_file_name_mk = qq{$temp_dir/output_mk.cram.seqchksum};
      my $cram_seqchksum_fifo_name_mk = qq{$temp_dir/output_mk.cram.seqchksum.fifo};
      my $cram_crai_file_name_mk = qq{$temp_dir/output_mk.cram.crai};
      my $cram_crai_fifo_name_mk = qq{$temp_dir/output_mk.cram.crai.fifo};
      my $cram_md5_file_name_mk = qq{$temp_dir/output_mk.cram.md5};
      my $cram_md5_fifo_name_mk = qq{$temp_dir/output_mk.cram.md5.fifo};
      my $cram_file_name_mk = qq{$temp_dir/output_mk.cram};
      my $cram_fifo_name_mk = qq{$temp_dir/output_mk.cram.fifo};
      my $bam_file_name_mk = qq{$temp_dir/output_mk.bam};
      my $bam_seqchksum_file_name_mk = qq{$temp_dir/output_mk.bam.seqchksum};
      my $bam_seqchksum_fifo_name_mk = qq{$temp_dir/output_mk.bam.seqchksum.fifo};

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-" }. qq{ > $temp_dir/output_mk.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_mk.flagstat};

      my $expected_bamcheck_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.bamcheck.fifo | };
      $expected_bamcheck_cmd .= qq{$bam_bamcheck_cmd > $temp_dir/output_mk.bamcheck};

      my $expected_index_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.index.fifo | /software/solexa/pkg/samtools/samtools-$samtools_version/samtools index /dev/stdin /dev/stdout > $temp_dir/output_mk.bai};

      my $expected_pb_cal_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.pb_cal.fifo | };
      $expected_pb_cal_cmd .= qq{$bam_pb_cal_cmd -p $temp_dir/output_mk -filter-bad-tiles 2 -};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ | tee $temp_dir/output_mk.bam.seqchksum.fifo > $temp_dir/output_mk.bam.seqchksum};
  
      my $expected_altchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.alt.bschk.fifo | };
      $expected_altchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam hash=sha512primesums512};
      $expected_altchksum_cmd .= qq{ > $temp_dir/output_mk.bam.sha512primesums512.seqchksum};
  
      my $expected_scramble_cmd = qq{$scramble -I bam -O cram -r t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta < $temp_dir/output_mk.bam.scramble.fifo };
      $expected_scramble_cmd .= qq{| tee $temp_dir/output_mk.cram.fifo $cram_crai_fifo_name_mk $cram_md5_fifo_name_mk > $temp_dir/output_mk.cram};

      my $expected_cramchksum_cmd =  qq{set -o pipefail; cat $cram_fifo_name_mk | $bamseqchksum verbose=0 inputformat=cram reference=t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta };
      $expected_cramchksum_cmd .= qq{| tee $cram_seqchksum_fifo_name_mk > $cram_seqchksum_file_name_mk};

      my $expected_cramindex_cmd = qq{set -o pipefail; cat $cram_crai_fifo_name_mk | $cram_index - $cram_crai_file_name_mk};
      my $expected_crammd5_cmd = qq{set -o pipefail; cat $cram_md5_fifo_name_mk | md5sum -b | tr -d }.q{"\n *-" }. qq{ > $cram_md5_file_name_mk};

      my $expected_diff_cmd = qq{diff $bam_seqchksum_fifo_name_mk $cram_seqchksum_fifo_name_mk};

      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_altchksum_cmd;
      push  @expected_fork_cmds, $expected_scramble_cmd;
      push  @expected_fork_cmds, $expected_cramchksum_cmd;
      push  @expected_fork_cmds, $expected_cramindex_cmd;
      push  @expected_fork_cmds, $expected_crammd5_cmd;
      push  @expected_fork_cmds, $expected_diff_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_bamcheck_cmd;
      push  @expected_fork_cmds, $expected_index_cmd;
      push  @expected_fork_cmds, $expected_pb_cal_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for fork generated correctly') or diag explain [$bam->fork_cmds(),$expected_fork_cmds];

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

      is (!-z "$temp_dir/output.bam", 1, 'BAM file created with contents');      
      is (!-z "$temp_dir/output.bai", 1, 'BAM index created with contents');      
      is (!-z "$temp_dir/output.bam.md5", 1, 'BAM md5 created with contents');      
      is (!-z "$temp_dir/output.flagstat", 1, 'BAM flagstat created with contents');      
      is (!-z "$temp_dir/output_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents');
      is (!-z "$temp_dir/output_quality_cycle_surv.txt", 1, 'Quality surv created with contents');
      is (!-z "$temp_dir/output_quality_error.txt", 1, 'Quality error table created with contents');
      is (!-z "$temp_dir/output.cram", 1, 'CRAM file created with contents');
      is (!-z "$temp_dir/output.cram.crai", 1, 'CRAM index file created with contents');
      is (!-z "$temp_dir/output.cram.md5", 1, 'CRAM md5 file created with contents');
      is (!-z "$temp_dir/output.seqchksum", 1, 'BAM seqchksum file created with contents');
      is (!-z "$temp_dir/output.sha512primesums512.seqchksum", 1, 'sha512primesums512 seqchksum file created with contents');
      is (!-e "$temp_dir/output_mk.cram.seqchksum", 1, 'CRAM seqchksum file created with contents has been removed');
  }    
}

{
### 15156_1#54.bam, is a subset using DownsampleSam.jar PROBABILITY=0.01, from study 3123 (current_studies.alignments_in_bam=0)
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
            40 unless ($ENV{'TOOLS_INSTALLED'});
      my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                {
                  input_bam     => 't/data/sequence/15156_1#54.bam',
                  output_bam    => "$temp_dir/non_aligned_output.bam",
                  metrics_json  => "$temp_dir/non_aligned_metrics.json",
                  temp_dir      => $temp_dir,
                  metrics_file  =>  $temp_dir . '/non_aligned_metrics.txt',
                  default_java_xmx_elc => $elc_memory_for_deployment,
                  default_java_xmx_bts => $bts_memory_for_deployment,
                });
      my $expected_mark_duplicate_cmd = qq{$bammarkduplicates I=t/data/sequence/15156_1#54.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/non_aligned_metrics.txt level=0};
      is($bam->mark_duplicate_cmd(), $expected_mark_duplicate_cmd, 'correct biobambam command');
      is ( $bam->no_alignment(), 1, 'input bam with no alignment');
      lives_ok{$bam->process()} q{Processed OK};                   
      is (-e "$temp_dir/non_aligned_output.cram", 1, 'non-aligned CRAM file created');
      is (-e "$temp_dir/non_aligned_output.cram.md5", 1, 'non-aligned CRAM md5 file created');
      is (-e "$temp_dir/non_aligned_output.cram.seqchksum", 1, 'non-aligned CRAM seqchksum file created');
      is (-e "$temp_dir/non_aligned_output.bam.seqchksum", 1, 'non-aligned BAM seqchksum file created');
      is (-e "$temp_dir/non_aligned_output.bam.sha512primesums512.seqchksum", 1, 'non-aligned sha512primesums512 seqchksum file created');
                            
      $bam = npg_common::sequence::BAM_MarkDuplicate->new(
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
      $expected_mark_duplicate_cmd = qq{$bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics_no_align.txt level=0};
      is($bam->mark_duplicate_cmd(), $expected_mark_duplicate_cmd, 'correct biobambam command');
      ok( $bam->no_alignment(), 'input bam without alignment');
      my $expected_bam_tag_stripper_cmd = qq{INPUT=t/data/sequence/unaligned.bam OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      my $bam_tag_stripper_cmd = $bam->bam_tag_stripper_cmd();
      like($bam_tag_stripper_cmd, qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      # stop qr// in like interpolating READ_NAME_REGEX by enclosing value in \Q..\E
      my $expected_elc_cmd = qq{INPUT=t/data/sequence/unaligned.bam OUTPUT=$temp_dir/metrics_no_align.txt TMP_DIR=$temp_dir READ_NAME_REGEX='\Q[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\E' VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'};
      like($bam->estimate_library_complexity_cmd(), qr/$expected_elc_cmd/, 'correct elc command');

      like($bam->bamseqchksum_cmd(q{bam}), qr{\Qbamseqchksum verbose=0 inputformat=bam\E}, 'correct bamseqchksum command for a bam file with reference but no alignment');
      like($bam->bamseqchksum_cmd(q{cram}), qr{\Qbamseqchksum verbose=0 inputformat=cram\E}, 'correct bamseqchksum command for a cram file with reference but no alignment');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;

      is ($bam->_result->info->{'Picard-tools'}, undef, 'test picard version is not defined if no_alignment flag used');

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file');

      my $expected_tee_cmd = qq{set -o pipefail;$bam_tag_stripper_cmd | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_no_align.bam.md5.fifo $temp_dir/output_no_align.bam.flagstat.fifo $temp_dir/output_no_align.bam.bamcheck.fifo $temp_dir/output_no_align.bam.bschk.fifo $temp_dir/output_no_align.bam.alt.bschk.fifo $temp_dir/output_no_align.bam.scramble.fifo > $temp_dir/output_no_align.bam};
      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly if no_alignment flag used');

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-" }. qq{ > $temp_dir/output_no_align.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_no_align.flagstat};

      my $expected_bamcheck_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.bamcheck.fifo | };
      $expected_bamcheck_cmd .= qq{$bam_bamcheck_cmd > $temp_dir/output_no_align.bamcheck};

      my $expected_altchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.alt.bschk.fifo | };
      $expected_altchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam hash=sha512primesums512};
      $expected_altchksum_cmd .= qq{ > $temp_dir/output_no_align.bam.sha512primesums512.seqchksum};

      my $cram_md5_file_name_mk = qq{$temp_dir/output_no_align.cram.md5};
      my $cram_md5_fifo_name_mk = qq{$temp_dir/output_no_align.cram.md5.fifo};
      my $cram_seqchksum_file_name_mk = qq{$temp_dir/output_no_align.cram.seqchksum};
      my $cram_seqchksum_fifo_name_mk = qq{$temp_dir/output_no_align.cram.seqchksum.fifo};
  
      my $expected_scramble_cmd = qq{$scramble -I bam -O cram < $temp_dir/output_no_align.bam.scramble.fifo };
      $expected_scramble_cmd .= qq{| tee $temp_dir/output_no_align.cram.fifo $cram_md5_fifo_name_mk > $temp_dir/output_no_align.cram};

      my $expected_cramchksum_cmd =  qq{set -o pipefail; cat $temp_dir/output_no_align.cram.fifo | $bamseqchksum verbose=0 inputformat=cram };
      $expected_cramchksum_cmd .= qq{| tee $cram_seqchksum_fifo_name_mk > $cram_seqchksum_file_name_mk};

      my $expected_crammd5_cmd = qq{set -o pipefail; cat $cram_md5_fifo_name_mk | md5sum -b | tr -d }.q{"\n *-" }. qq{ > $cram_md5_file_name_mk};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ | tee $temp_dir/output_no_align.bam.seqchksum.fifo};
      $expected_bamchksum_cmd .= qq{ > $temp_dir/output_no_align.bam.seqchksum};
  
      my $expected_diff_cmd = qq{diff $temp_dir/output_no_align.bam.seqchksum.fifo $cram_seqchksum_fifo_name_mk};

      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_altchksum_cmd;
      push  @expected_fork_cmds, $expected_scramble_cmd;
      push  @expected_fork_cmds, $expected_cramchksum_cmd;
      push  @expected_fork_cmds, $expected_crammd5_cmd;
      push  @expected_fork_cmds, $expected_diff_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_bamcheck_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for ForkManager generated correctly') or diag explain [$bam->fork_cmds(),$expected_fork_cmds];

      lives_ok{$bam->process()} q{Processed OK};
      is (-e "$temp_dir/output_no_align.bam.md5.fifo", 1, 'md5 FIFO created');
      is (-e "$temp_dir/output_no_align.bam.flagstat.fifo", 1, 'flagstat FIFO created');
      is (-e "$temp_dir/output_no_align.bam.bamcheck.fifo", 1, 'bamcheck FIFO created');
      is (!-e "$temp_dir/output_no_align.bam.index.fifo", 1, 'index FIFO NOT created');
      is (!-e "$temp_dir/output_no_align.bam.pb_cal.fifo", 1, 'pb_cal FIFO NOT created');
      is (-e "$temp_dir/output_no_align.bam.scramble.fifo", 1, 'scramble FIFO created');
      is (-e "$temp_dir/output_no_align.bam.bschk.fifo", 1, 'bamseqchksum input FIFO created');
      is (-e "$temp_dir/output_no_align.bam.seqchksum.fifo", 1, 'bamseqchksum output FIFO created');
      is (-e "$temp_dir/output_no_align.cram.seqchksum.fifo", 1, 'bamseqchksum output FIFO created');

      is (!-z "$temp_dir/output_no_align.bam", 1, 'BAM file created with contents');      
      is (!-e "$temp_dir/output_no_align.bai", 1, 'BAM index NOT created if no_alignment flag used');      
      is (!-z "$temp_dir/output_no_align.bam.md5", 1, 'BAM md5 created with contents');      
      is (!-z "$temp_dir/output_no_align.flagstat", 1, 'BAM flagstat created with contents');      
      is (!-z "$temp_dir/output_no_align_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents');
      is (!-z "$temp_dir/output_no_align_quality_cycle_surv.txt", 1, 'Quality surv created with contents');
      is (!-z "$temp_dir/output_no_align_quality_error.txt", 1, 'Quality error table created with contents');
      is (-e "$temp_dir/output_no_align.cram", 1, 'CRAM file created if no_alignment flag used');
      is (-e "$temp_dir/output_no_align.cram.seqchksum", 1, 'CRAM seqchksum file created if no_alignment flag used');
      is (!-z "$temp_dir/output_no_align.bam.seqchksum", 1, 'BAM seqchksum file created with contents if no_alignment flag used');
  }
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         32 unless ($ENV{'TOOLS_INSTALLED'});
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
                 replace_file => 0,
               });
      my $expected_mark_duplicate_cmd = qq{$bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics_phix.txt level=0};
      is($bam->mark_duplicate_cmd(), $expected_mark_duplicate_cmd, 'correct biobambam command');
      ok(!$bam->no_alignment(), 'input PhiX bam with alignment');
      
      my $expected_bam_tag_stripper_cmd = qq{INPUT=/dev/stdin OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      my $bam_tag_stripper_cmd = $bam->bam_tag_stripper_cmd();
      like($bam_tag_stripper_cmd, qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      like($bam->bamseqchksum_cmd(q{bam}), qr{\Qbamseqchksum verbose=0 inputformat=bam\E}, 'correct bamseqchksum command for a bam file with reference but no alignment');
      like($bam->bamseqchksum_cmd(q{cram}), qr{\Qbamseqchksum verbose=0 inputformat=cram\E}, 'correct bamseqchksum command for a cram file with reference but no alignment');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file');

      my $expected_tee_cmd = qq{set -o pipefail;$bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics_phix.txt level=0 | $bam_tag_stripper_cmd | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_phix.bam.md5.fifo $temp_dir/output_phix.bam.flagstat.fifo $temp_dir/output_phix.bam.bamcheck.fifo $temp_dir/output_phix.bam.bschk.fifo $temp_dir/output_phix.bam.alt.bschk.fifo $temp_dir/output_phix.bam.index.fifo $temp_dir/output_phix.bam.pb_cal.fifo $temp_dir/output_phix.bam.scramble.fifo > $temp_dir/output_phix.bam};
      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly for PhiX');

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-" }. qq{ > $temp_dir/output_phix.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_phix.flagstat};

      my $expected_bamcheck_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.bamcheck.fifo | };
      $expected_bamcheck_cmd .= qq{$bam_bamcheck_cmd > $temp_dir/output_phix.bamcheck};

      my $expected_index_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.index.fifo | /software/solexa/pkg/samtools/samtools-$samtools_version/samtools index /dev/stdin /dev/stdout > $temp_dir/output_phix.bai};

      my $bam_pb_cal_cmd = $bam->pb_cal_cmd();
      my $expected_pb_cal_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.pb_cal.fifo | };
      $expected_pb_cal_cmd .= qq{$bam_pb_cal_cmd -p $temp_dir/output_phix -filter-bad-tiles 2 -};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ | tee $temp_dir/output_phix.bam.seqchksum.fifo};
      $expected_bamchksum_cmd .= qq{ > $temp_dir/output_phix.bam.seqchksum};
 
      my $cram_seqchksum_file_name_phix = qq{$temp_dir/output_phix.cram.seqchksum};
      my $cram_seqchksum_fifo_name_phix = qq{$temp_dir/output_phix.cram.seqchksum.fifo};
      my $cram_crai_file_name_phix = qq{$temp_dir/output_phix.cram.crai};
      my $cram_crai_fifo_name_phix = qq{$temp_dir/output_phix.cram.crai.fifo};
      my $cram_md5_file_name_phix = qq{$temp_dir/output_phix.cram.md5};
      my $cram_md5_fifo_name_phix = qq{$temp_dir/output_phix.cram.md5.fifo};
      my $cram_file_name_phix = qq{$temp_dir/output_phix.cram};
      my $cram_fifo_name_phix = qq{$temp_dir/output_phix.cram.fifo};

      my $expected_altchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.alt.bschk.fifo | };
      $expected_altchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam hash=sha512primesums512};
      $expected_altchksum_cmd .= qq{ > $temp_dir/output_phix.bam.sha512primesums512.seqchksum};
  
      my $expected_scramble_cmd = qq{$scramble -I bam -O cram < $temp_dir/output_phix.bam.scramble.fifo };
      $expected_scramble_cmd .= qq{| tee $temp_dir/output_phix.cram.fifo $cram_crai_fifo_name_phix $cram_md5_fifo_name_phix > $temp_dir/output_phix.cram};

      my $expected_cramchksum_cmd =  qq{set -o pipefail; cat $cram_fifo_name_phix | $bamseqchksum verbose=0 inputformat=cram };
      $expected_cramchksum_cmd .= qq{| tee $cram_seqchksum_fifo_name_phix > $cram_seqchksum_file_name_phix};

      my $expected_cramindex_cmd = qq{set -o pipefail; cat $cram_crai_fifo_name_phix | $cram_index - $cram_crai_file_name_phix};
      my $expected_crammd5_cmd = qq{set -o pipefail; cat $cram_md5_fifo_name_phix | md5sum -b | tr -d }.q{"\n *-" }. qq{ > $cram_md5_file_name_phix};

      my $expected_diff_cmd = qq{diff $temp_dir/output_phix.bam.seqchksum.fifo $cram_seqchksum_fifo_name_phix};

      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_altchksum_cmd;
      push  @expected_fork_cmds, $expected_scramble_cmd;
      push  @expected_fork_cmds, $expected_cramchksum_cmd;
      push  @expected_fork_cmds, $expected_cramindex_cmd;
      push  @expected_fork_cmds, $expected_crammd5_cmd;
      push  @expected_fork_cmds, $expected_diff_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_bamcheck_cmd;
      push  @expected_fork_cmds, $expected_index_cmd;
      push  @expected_fork_cmds, $expected_pb_cal_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for ForkManager generated correctly') or diag explain [$bam->fork_cmds(),$expected_fork_cmds];

      lives_ok{$bam->process()} q{Processed OK};

      is (-e "$temp_dir/output_phix.bam.md5.fifo", 1, 'md5 FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.flagstat.fifo", 1, 'flagstat FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.bamcheck.fifo", 1, 'bamcheck FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.index.fifo", 1, 'index FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.pb_cal.fifo", 1, 'pb_cal FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.scramble.fifo", 1, 'scramble FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.bschk.fifo", 1, 'bamseqchksum input FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.bam.seqchksum.fifo", 1, 'bamseqchksum output FIFO created for PhiX');
      is (-e "$temp_dir/output_phix.cram.seqchksum.fifo", 1, 'bamseqchksum output FIFO created for PhiX');

      is (!-z "$temp_dir/output_phix.bam", 1, 'BAM file created with contents for PhiX');      
      is (!-z "$temp_dir/output_phix.bai", 1, 'BAM index created with contents for PhiX');      
      is (!-z "$temp_dir/metrics_phix.bam.json", 1, 'metrics json created with contents for PhiX');      
      is (!-z "$temp_dir/metrics_phix.txt", 1, 'metrics txt created with contents for PhiX');      
      is (!-z "$temp_dir/output_phix.bam.md5", 1, 'BAM md5 created with contents for PhiX');      
      is (!-z "$temp_dir/output_phix.flagstat", 1, 'BAM flagstat created with contents for PhiX');      
      is (!-z "$temp_dir/13388_2#40_phix_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents for PhiX');
      is (!-z "$temp_dir/13388_2#40_phix_quality_cycle_surv.txt", 1, 'Quality surv created with contents for PhiX');
      is (!-z "$temp_dir/13388_2#40_phix_quality_error.txt", 1, 'Quality error table created with contents for PhiX');

      is (!-z "$temp_dir/output_no_align.cram", 1, 'CRAM file created with contents for PhiX');
      is (!-z "$temp_dir/output_no_align.cram.seqchksum", 1, 'CRAM seqchksum file created with contents for PhiX');
      is (!-z "$temp_dir/output_no_align.bam.seqchksum", 1, 'BAM seqchksum file created with contents for PhiX');
  }
}

1;
