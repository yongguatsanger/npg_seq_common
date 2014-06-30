#########
# Author:        gq1
# Created:       2009-06-21
#

use strict;
use warnings;
use Test::More tests => 46;
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
            9 unless ($ENV{'TOOLS_INSTALLED'});
  my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam     => 'input.bam',
                 output_bam    => 'output.bam',
                 metrics_json  => 'metrics.json',
                 not_strip_bam_tag => 1,
                 no_alignment => 0,
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
  is($bam->bamseqchksum_cmd(q{bam}), q{bamseqchksum verbose=1 inputformat=bam > output.bam.seqchksum}, 'correct bamseqchsum command for a bam file');
  is($bam->bamseqchksum_cmd(q{cram}), q{bamseqchksum verbose=1 inputformat=cram > output.cram.seqchksum}, 'correct bamseqchsum command for a cram file with no reference');
  
       };
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         19 unless ($ENV{'TOOLS_INSTALLED'});
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
               {
                 input_bam     => 't/data/sequence/SecondCall/4392_1.bam',
                 output_bam    => "$temp_dir/output_mk.bam",
                 metrics_json  => "$temp_dir/metrics.json",
                 sort_input    => 1,
                 temp_dir      => $temp_dir,
                 metrics_file  => $temp_dir . '/metrics.txt',
                 not_strip_bam_tag => 1,
                 reference => 't/data/references/Plasmodium_falciparum/3D7_Oct11v3/all/fasta/Pf3D7_v3.fasta',
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
      is($bam_pb_cal_cmd, q{/software/solexa/pkg/pb_calibration/v10.14/bin/calibration_pu}, 'correct pb_cal comand');

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file, using newer samtools than current');

      my $expected_bamseqchk_cmd = qq{bamseqchksum verbose=1 inputformat=cram reference=t/data/references/Plasmodium_falciparum/3D7_Oct11v3/all/fasta/Pf3D7_v3.fasta > $temp_dir/output_mk.cram.seqchksum};
      is($bam->bamseqchksum_cmd(q{cram}), $expected_bamseqchk_cmd, 'correct bamseqchsum command for a cram file with reference');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      ok ($bam->_result->info->{'Samtools'}, 'samtools version is defined for an unaligned bam');
      ok ($bam->_result->info->{'Picard-tools'}, 'test picard version is defined for an unaligned bam');

      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;
      lives_ok{$bam->process()} q{Processed OK};

      my $expected_tee_cmd = qq{set -o pipefail;/software/hpag/biobambam/0.0.147/bin/bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics.txt level='0' | $bam_tag_stripper_cmd | tee  >(md5sum -b | tr -d }.q{"\n *-"};
      $expected_tee_cmd .= qq{ > $temp_dir/output_mk.bam.md5) >(/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_mk.flagstat) >($bam_bamcheck_cmd > $temp_dir/output_mk.bamcheck) >(/software/solexa/pkg/samtools/samtools-$samtools_version/samtools index /dev/stdin /dev/stdout > $temp_dir/output_mk.bai) };
      $expected_tee_cmd .= qq{>(/software/badger/bin/scramble -I bam -O cram -r t/data/references/Plasmodium_falciparum/3D7_Oct11v3/all/fasta/Pf3D7_v3.fasta  | tee >(bamseqchksum verbose=1 inputformat=cram reference=t/data/references/Plasmodium_falciparum/3D7_Oct11v3/all/fasta/Pf3D7_v3.fasta > $temp_dir/output_mk.cram.seqchksum) > $temp_dir/output_mk.cram ) };
      $expected_tee_cmd .= qq{>($bam_pb_cal_cmd -p t/data/sequence/plasmodium -filter-bad-tiles 2 -) >(bamseqchksum verbose=1 inputformat=bam > $temp_dir/output_mk.bam.seqchksum)  > $temp_dir/output_mk.bam};

      is($bam->tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly');
      is (-e "$temp_dir/output_mk.bam", 1, 'BAM file created');      
      is (-e "$temp_dir/output_mk.bai", 1, 'BAM index created');      
      is (-e "$temp_dir/output_mk.bam.md5", 1, 'BAM md5 created');      
      is (-e "$temp_dir/output_mk.flagstat", 1, 'BAM flagstat created');      
      is (-e "$temp_dir/output_mk.cram", 1, 'CRAM file created');
  }    
}

{
  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         17 unless ($ENV{'TOOLS_INSTALLED'});
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
               {
                 input_bam     => 't/data/sequence/unaligned.bam',
                 output_bam    => "$temp_dir/output_no_align.bam",
                 metrics_json  => "$temp_dir/metrics.json",
                 sort_input    => 1,
                 temp_dir      => $temp_dir,
                 metrics_file  => $temp_dir . '/metrics.txt',
                 not_strip_bam_tag => 0,
                 default_java_xmx_elc => $elc_memory_for_deployment,
                 default_java_xmx_bts => $bts_memory_for_deployment,
               });
      my $expected_mark_duplicate_cmd = qq{bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$temp_dir/metrics.txt};
      like($bam->mark_duplicate_cmd(), qr/$expected_mark_duplicate_cmd/, 'correct biobambam command');
      ok( $bam->no_alignment(), 'input bam without alignment');
      my $expected_bam_tag_stripper_cmd = qq{INPUT=t/data/sequence/unaligned.bam OUTPUT=/dev/stdout TMP_DIR=$temp_dir CREATE_INDEX='FALSE' CREATE_MD5_FILE='FALSE' VALIDATION_STRINGENCY='SILENT' VERBOSITY='INFO' STRIP='OQ' KEEP='a3' KEEP='aa' KEEP='af' KEEP='ah' KEEP='as' KEEP='br' KEEP='qr' KEEP='tq' KEEP='tr'};
      my $bam_tag_stripper_cmd = $bam->bam_tag_stripper_cmd();
      like($bam_tag_stripper_cmd, qr/$expected_bam_tag_stripper_cmd/, 'correct bam_tag_stripper command');
      
      # stop qr// in like interpolating READ_NAME_REGEX by enclosing value in \Q..\E
      my $expected_elc_cmd = qq{INPUT=t/data/sequence/unaligned.bam OUTPUT=$temp_dir/metrics.txt TMP_DIR=$temp_dir READ_NAME_REGEX='\Q[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\E' VALIDATION_STRINGENCY='SILENT' VERBOSITY='ERROR'};
      like($bam->estimate_library_complexity_cmd(), qr/$expected_elc_cmd/, 'correct elc command');

      is($bam->bamseqchksum_cmd(q{bam}), qq{bamseqchksum verbose=1 inputformat=bam > $temp_dir/output_no_align.bam.seqchksum}, 'correct bamseqchsum command for a bam file with reference but no alignment');
      is($bam->bamseqchksum_cmd(q{cram}), qq{bamseqchksum verbose=1 inputformat=cram > $temp_dir/output_no_align.cram.seqchksum}, 'correct bamseqchsum command for a cram file with reference but no alignment');

      lives_ok {$bam->_version_info} 'getting tools version info lives';
      my $samtools_version_str = $bam->_result->info->{'Samtools'};
      ok ($samtools_version_str, 'samtools version is defined');
      my ($samtools_version, $samtools_revison) = split / /, $samtools_version_str;

      is ($bam->_result->info->{'Picard-tools'}, undef, 'test picard version is not defined if no_alignment flag used');

      my $bam_bamcheck_cmd = $bam->bamcheck_cmd();
      is($bam_bamcheck_cmd, q{/software/solexa/pkg/samtools/samtools-0.1.19/misc/bamcheck}, 'correct bamcheck command for bam file');

      lives_ok{$bam->process()} q{Processed OK};

      my $expected_tee_cmd = qq{set -o pipefail;$bam_tag_stripper_cmd | tee  >(md5sum -b | tr -d }.q{"\n *-"}. qq{ > $temp_dir/output_no_align.bam.md5) >(/software/solexa/pkg/samtools/samtools-$samtools_version/samtools flagstat -  > $temp_dir/output_no_align.flagstat) >($bam_bamcheck_cmd > $temp_dir/output_no_align.bamcheck) >(bamseqchksum verbose=1 inputformat=bam > $temp_dir/output_no_align.bam.seqchksum)  > $temp_dir/output_no_align.bam};
      is($bam->tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly if no_alignment flag used');
      is (-e "$temp_dir/output_no_align.bam", 1, 'BAM file created');      
      is (!-e "$temp_dir/output_no_align.bai", 1, 'BAM index NOT created if no_alignment flag used');      
      is (-e "$temp_dir/output_no_align.bam.md5", 1, 'BAM md5 created');      
      is (-e "$temp_dir/output_no_align.flagstat", 1, 'BAM flagstat created');      
      is (!-e "$temp_dir/output_no_align.cram", 1, 'CRAM file NOT created if no_alignment flag used');
  }    
}


1;
