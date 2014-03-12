#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2011-08-02
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$

package npg_common::bam_align_irods;

use Moose;
use Carp;
use English qw(-no_match_vars);
use File::Basename;
use File::Copy;
use File::Spec::Functions qw(catfile);
use File::Temp qw( tempfile tempdir );
use Perl6::Slurp;
use IPC::Open3;
use Cwd qw(abs_path);
use JSON;
use FindBin qw($Bin);

use npg_common::irods::Loader;
use npg_common::irods::run::Bam;
use npg_tracking::data::reference;
use npg_qc::Schema;

with qw(
    MooseX::Getopt
    npg_common::roles::log
    npg_common::roles::software_location
);

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/msx; $r; };

Readonly::Scalar my $BAM_ALIGNMENT_CMD         => 'bam_alignment.pl';
Readonly::Scalar my $BAM_FLAGSTATS_LOADER_CMD  => 'npg_qc_autoqc_data.pl';

Readonly::Scalar my $CHANGE_BAM_HEADER_JAR     => 'ChangeBamHeader.jar';
Readonly::Scalar my $STRIP_BAM_TAG_JAR         => 'BamTagStripper.jar';

Readonly::Scalar my $DEFAULT_TAG_INDEX => -1;
Readonly::Scalar our $EXIT_CODE_SHIFT => 8;

has 'change_bam_header_jar' => (
              is            => 'ro',
              isa           => 'NpgCommonResolvedPathJarFile',
              lazy_build    => 1,
              coerce        => 1,
              documentation => 'full path to BamHeader.jar',
    );
sub _build_change_bam_header_jar { return $CHANGE_BAM_HEADER_JAR; }

has 'strip_bam_tag_jar'     => (
              is            => 'ro',
              isa           => 'NpgCommonResolvedPathJarFile',
              lazy_build    => 1,
              coerce        => 1,
              documentation => 'full path to StripBamTag.jar',
    );
sub _build_strip_bam_tag_jar { return $STRIP_BAM_TAG_JAR; }

has 'is_paired_read' => (
              is            => 'ro',
              isa           => 'Bool',
              default       => 1,
              documentation => 'is_paired_read flag, true by default, set to false manually if needed',
    );

has 'spiked_phix_split' => (
              is            => 'ro',
              isa           => 'Bool',
              default       => 0,
              documentation => 'spiked_phix_split flag, false by default, set to true manually if needed',
    );

has 'non_consent_split' => (
              is            => 'ro',
              isa           => 'Bool',
              default       => 0,
              documentation => 'non_consent_split flag, false by default, set to true manually if needed',
    );

has 'contains_nonconsented_xahuman' => (
              is            => 'ro',
              isa           => 'Bool',
              default       => 0,
              documentation => 'contains_nonconsented_xahuman flag, false by default, set to true manually if needed',
    );

has 'separate_y_chromosome_data' => (
              is            => 'ro',
              isa           => 'Bool',
              default       => 0,
              documentation => 'separate_y_chromosome_data flag, false by default, set to true manually if needed',
    );

has 'bamsort_cmd'   => ( is      => 'ro',
                         isa     => 'NpgCommonResolvedPathExecutable',
                         coerce  => 1,
                         default => 'bamsort',
                       );

has 'bamcollate2_path'  => (
                       is      => 'ro',
                       isa     => 'NpgCommonResolvedPathExecutable',
                       lazy    => 1,
                       builder => '_build_bamcollate2_path',
                              );
sub _build_bamcollate2_path {
  my $self = shift;
  return catfile($self->biobambam_bin,'bamcollate2');
}

has 'biobambam_bin' => (
                       isa     => 'NpgTrackingDirectory',
                       is      => 'ro',
                       lazy    => 1,
                       builder => '_build_biobambam_bin',
                            );
sub _build_biobambam_bin {
  my $self = shift;
  return dirname($self->bamsort_cmd());
}


sub run {
   my $self = shift;
   $self->qc_dir();
   my $task = $self->task();

   $self->get_old_bam_from_irods();

   if( $task eq q{bam_realign} ){
      $self->collate_old_bam();
   }

   $self->$task();

   $self->reload_bam();

   if( $task eq q{bam_realign} ){
      my $count = 0;

      if( $self->spiked_phix_split ){
         $self->reload_bam('phix');
         $count++;
      }

      if( $self->non_consent_split ){
         $self->reload_bam('human');
         $count++;
      }

      if( $self->contains_nonconsented_xahuman ){
         $self->reload_bam('xahuman');
         $count++;
      }

      if( $self->separate_y_chromosome_data ){
         $count++;
         $self->reload_bam('yhuman');
      }

      if( $count ) {
         $self->convert_alignment_filter_json();
      }

      # delete pre-existing flagstats for the input bam file only
      $self->delete_qc_bam_flagstats();

      # if the bam_realign produced multiple flagstats files we only reload those
      # which were deleted by delete_qc_bam_flagstats() or were already in the db
      $self->reload_bam_flagstats();
   }

   return;
}

has 'id_run' => (
    is            => 'ro',
    isa           => 'Int',
    required      => 1,
);

has 'lane' => (
    is            => 'ro',
    isa           => 'Int',
    required      => 1,
);

has 'tag_index' => (
    is            => 'ro',
    isa           => 'Int',
    required      => 0,
);

has 'alignment_filter' => (isa        => 'Str',
                           is         => 'rw',
                           required   => 0,
                          );

has 'rt_ticket' => (
    is            => 'ro',
    isa           => 'Int',
    required      => 0,
);

has 'task' => (
    is            => 'ro',
    isa           => 'Str',
    required      => 1,
);

has 'reference' => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
);

sub _build_reference {
   my $self = shift;
   my $ref_finder = npg_tracking::data::reference->new(
                                        id_run   => $self->id_run(),
                                        position => $self->lane(),
                                        tag_index => $self->tag_index(),
   );
   my $refs = $ref_finder->refs;
   if(scalar @{$refs} >1){
     croak 'more than one reference found';
   }
   return $refs->[0];
}

has 'human_reference'   => (isa           => 'Str',
                            is            => 'rw',
                            required      => 0,
                            lazy_build    => 1,
                            documentation => 'A path to human reference for bwa alignment to split non-consent reads',
                           );

sub _build_human_reference {
   my $self = shift;

   my $ruser = Moose::Meta::Class->create_anon_class(
          roles => [qw/npg_tracking::data::reference::find/])
          ->new_object({species => q{Homo_sapiens}});

   my $human_ref;
   eval {
      $human_ref = $ruser->refs->[0];
      1;
   } or do{
      carp $EVAL_ERROR;
   };

   return $human_ref;
}

has  '_lims'    =>  ( isa             => 'Maybe[st::api::lims]',
                      is              => 'ro',
                      init_arg        => undef,
                      lazy_build      => 1,
                    );
sub _build__lims {
    my $self = shift;
    if( defined $self->id_run() ){
       return st::api::lims->new(
                            id_run  => $self->id_run(),
                            position  => $self->lane(),
                            tag_index => $self->tag_index(),
                                );
    }
    $self->log( 'id_run, position or tag_index need to be set to get lims information' );
    return;
}

has 'working_dir' => (
    is            => 'ro',
    isa           => 'Str',
    required      => 1,
);

has 'bam_filename' => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
);
sub _build_bam_filename {
   my $self = shift;

   my $name = $self->id_run().q{_}.$self->lane();

   if(defined $self->tag_index() ){
      $name .= q{#}.$self->tag_index();
   }
   if( $self->alignment_filter() ){
      $name .= q{_}.$self->alignment_filter();
   }
   return $name.q{.bam};
}

has 'old_bam_filename' => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
);
sub _build_old_bam_filename {
   my $self = shift;

   my $bam_filename = $self->bam_filename();

   $bam_filename =~ s/[.]bam$/_old\.bam/mxs;

   return $bam_filename;
}

has 'qc_dir' => (
    is            => 'ro',
    isa           => 'Str',
    lazy_build    => 1,
);
sub _build_qc_dir {
   my $self = shift;

   my $qc_dir = $self->working_dir().q{/qc/};
   system "mkdir -p $qc_dir";
   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
       croak "create qc directory: $ERRNO\nqc_dir";
   }
   return $qc_dir;
}

sub get_old_bam_from_irods {
   my $self = shift;

   mkdir $self->working_dir();
   my $command = q{iget -f -K /seq/}.$self->id_run().q{/}.$self->bam_filename(). q{ }. $self->working_dir().q{/}.$self->old_bam_filename();
   $self->log($command);
   my $iget_rs = system $command;
   if($iget_rs){
      croak "iget failed: $command";
   }
   return;
}

sub collate_old_bam {
   my $self = shift;

   my $old_bam = $self->working_dir().q{/}.$self->old_bam_filename();
   my $uncollated_old_bam = $self->working_dir().q{/}.$self->old_bam_filename();
   $uncollated_old_bam =~ s/[.]bam$/_uncollated.bam/mxs;

   $self->log("move $old_bam $uncollated_old_bam");
   move $old_bam, $uncollated_old_bam;

   # tags to keep
   # RG read group
   # PG program
   # BC barcode sequence
   # RT Sanger used to use RT for the barcode sequence
   # QT barcode quality
   # a3 adapter length
   # ah adapter match flag
   # tr transposon sequence
   # tq transposon quality
   # br inline index sequence
   # qr inline index quality

   my $collate_cmd = $self->biobambam_bin() . q{/} . q{bamsort SO=queryname level=0};
   $collate_cmd .= q{ tmpfile=} . tempdir( DIR => $self->working_dir(), CLEANUP => 1 ) . q{/} . q{bamsort};
   $collate_cmd .= q{ I=} . $uncollated_old_bam;
   $collate_cmd .= q{ | } . $self->biobambam_bin() . q{/} . q{bamcollate2 reset=1 resetaux=0 auxfilter=RG,PG,BC,RT,QT,a3,ah,tr,tq,br,qr};
   $collate_cmd .= q{ T=} . tempdir( DIR => $self->working_dir(), CLEANUP => 1 ) . q{/} . q{bamcollate2};
   $collate_cmd .= q{ > } . $old_bam;

   $self->log($collate_cmd);
   my $collate_cmd_rs = system $collate_cmd;
   if($collate_cmd_rs){
      croak "collate failed: $collate_cmd";
   }
   return;
}

sub bam_realign {
   my $self = shift;

   if( $self->alignment_filter && $self->alignment_filter ne 'nonhuman' ){
      croak "You cannot realign $self->alignment_filter bam file\n";
   }

   my $old_bam = $self->working_dir().q{/}.$self->old_bam_filename();
   my $output_prefix =  $self->working_dir().q{/}. $self->id_run().q{_}.$self->lane();
   if(defined $self->tag_index() ){
      $output_prefix .= q{#}.$self->tag_index();
   }
   my $reference = $self->reference();
   my $ref_dict = $reference;
   $ref_dict =~ s/bwa/picard/mxs;
   $ref_dict .= q{.dict};

   FindBin::again();
   my $command = catfile ($Bin, $BAM_ALIGNMENT_CMD)
               . qq{ --input $old_bam --output_prefix $output_prefix --reference $reference --ref_dict $ref_dict}
               . q{ --id_run } .$self->id_run() . q{ --position } . $self->lane();
   if( defined $self->tag_index() ){
       $command .= q{ --tag_index } . $self->tag_index();
   }
   if( $self->is_paired_read ) {
       $command  .= q{ --is_paired_read};
   }

   $command .= q{ --do_markduplicates };
   $command .= ($self->spiked_phix_split             ? q{ --} : q{ --no}) . q{spiked_phix_split};
   $command .= ($self->non_consent_split             ? q{ --} : q{ --no}) . q{non_consent_split};
   $command .= ($self->contains_nonconsented_xahuman ? q{ --} : q{ --no}) . q{contains_nonconsented_xahuman};
   $command .= ($self->separate_y_chromosome_data    ? q{ --} : q{ --no}) . q{separate_y_chromosome_data};

   $self->log($command);

   my $iget_rs = system "set -o pipefail; $command";
   if($iget_rs){
      croak "realignment failed: $command";
   }
   return;
}

sub reload_bam {
    my ( $self, $human_split ) = @_;


    my $bam_filename = $self->bam_filename();
    if( $human_split ) {
       $bam_filename =~ s/[.]bam$/_${human_split}.bam/mxs;
    }

    my $file = $self->working_dir().q{/}.$bam_filename;

    my $meta = {};
    if ( $human_split ){
       $meta->{target} = 0;
       if( $human_split ne 'phix' ) {
          $meta->{reference} = $self->human_reference();
       }
       if( $human_split eq 'yhuman' ) {
          $meta->{target} = 1;
       }
       $meta->{id_run} = $self->id_run();
       $meta->{lane}   = $self->lane();
       $meta->{type}   = q{bam};
       $meta->{alignment_filter} = $human_split;
       $meta->{alignment} = 1;

       my $total_reads = $self->_get_number_of_reads_qc($file);
       $meta->{total_reads} = $total_reads;

       # empty bam files are valid e.g. if no phix was spiked into the lane, but we should not
       # load an empty bam file better to re-run the script but not create the empty bam file
       my %human_split_options = {
                                  'phix'    => 'spiked_phix_split',
                                  'human'   => 'non_consent_split',
                                  'xahuman' => 'contains_nonconsented_xahuman',
                                  'yhuman'  => 'separate_y_chromosome_data',
                                 };
       croak "No reads in $file, re-run with --no".$human_split_options{$human_split} if $total_reads == 0;
    } else {
       $meta = {
                reference  => $self->reference(),
                alignment  => 1,
               };
    }

    if($self->rt_ticket()){
        $meta->{rt_ticket} = $self->rt_ticket();
    }

    my $loader = npg_common::irods::Loader->new({
       file         => $file,
       meta_data    => $meta,
       new_filename => $bam_filename,
       collection   => '/seq/'.$self->id_run(),
    });
    $loader->run();
    my $irods_file = $loader->collection().q{/}.$loader->new_filename();
    $self->remove_obsolete_replication($irods_file);

    my $bam_header = npg_common::irods::run::Bam->new(
                                          id_run => $self->id_run(),
                                          $self->resolved_paths,
                                      )->get_bam_reference($file);
    my $alignment = $bam_header->{alignment};
    if(!$alignment){
       $self->log('Bam file not aligned, no need to load bai file');
       return 1;
    }

    $meta = {type => 'bai'};
    $file =~ s/bam$/bai/mxs;
    my $irods_bai = $bam_filename;
    $irods_bai =~ s/bam$/bai/mxs;
    $loader = npg_common::irods::Loader->new({
       file         => $file,
       meta_data    => $meta,
       new_filename => $irods_bai,
       collection   => '/seq/'.$self->id_run(),
    });
    $loader->run();
    $irods_file = $loader->collection().q{/}.$loader->new_filename();
    $self->remove_obsolete_replication($irods_file);
    return;
}

sub remove_obsolete_replication {
   my ($self, $bam_file) = @_;

   my $rep_num_to_remove = $self->_get_rep_num_to_remove($bam_file);
   foreach my $rep_num (@{$rep_num_to_remove}){

       my $irm_command = qq{irm -n $rep_num $bam_file}; # irm on irods_path
       $self->log($irm_command);
       my $irm_rs =  system $irm_command;
       if($irm_rs){
          croak "Failed: $irm_command";
       }
   }
   return;
}

sub _get_rep_num_to_remove {
    my ($self, $bam_file) = @_;

    my @rep_num_to_remove = ();

    my $ils_cmd = qq{ils -l $bam_file}; # ils on irods_path

    my $pid = open3( undef, my $ils_out_fh, undef, $ils_cmd);

    while (my $line = <$ils_out_fh> ) {

       chomp $line;

       if ( $line =~ /ERROR[:][ ]lsUtil[:][ ]/mxs ){
          #file not exists in the server
          last;
       }elsif ( $line !~ /[ ]&[ ]/mxs ){
          my ($rep_num) = $line =~ /^[ ]{2}\w+\s+(\d+)/mxs;
          if(defined $rep_num){
             push @rep_num_to_remove, $rep_num;
          }
       }
    }

    waitpid $pid, 0;

    if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
       carp "Problems to list bam file $bam_file";
     }

    close $ils_out_fh or croak "can not close ils command output: $ERRNO";

    return \@rep_num_to_remove;
}


sub delete_qc_bam_flagstats {
   my $self = shift;

   my $schema = npg_qc::Schema->connect( );

   my $query = {id_run => $self->id_run(), position => $self->lane()};
   if( defined $self->tag_index() ){
       $query->{tag_index} = $self->tag_index();
   }else {
       $query->{tag_index} = $DEFAULT_TAG_INDEX;
   }
   if( $self->alignment_filter() ){
      $query->{human_split} = $self->alignment_filter();
   }else {
      $query->{human_split} = q{all};
   }
   $self->log('Deleting npgqc bam flagstats for run '
             . $query->{id_run}
             . q{ lane } . $query->{position}
             . q{ tag } . $query->{tag_index}
             . q{ human_split } . $query->{human_split}
           );
   $schema->txn_do(
       sub {
           my $result = $schema->resultset('BamFlagstats')
                         ->search( $query )->delete_all();
           if($result){
              $self->log('deletion successfull');
           }
       }
   );

   return;
}

sub reload_bam_flagstats {
   my $self = shift;
   FindBin::again();
   my $command = $BAM_FLAGSTATS_LOADER_CMD
               . q{ --id_run } . $self->id_run()
               . q{ --path } . $self->qc_dir();
   $self->log($command);
   my $loader_rs =  system $command;
   if($loader_rs){
      croak "reload bam flagstats failed: $command";
   }
   return;
}

sub change_bam_header {
    my $self = shift;

    $self->log('geting sample, library, and study informtation');

    my $lims = $self->_lims();
    my $sample_name = $lims->sample_publishable_name();
    my $library_id = $lims->library_id();
    my $study_name = $lims->study_publishable_name();
    my $study_description = $lims->study_description();

    my $old_bam = $self->working_dir().q{/}.$self->old_bam_filename();
    my $new_bam = $self->working_dir().q{/}.$self->bam_filename();

    my $change_header_cmd = $self->java_cmd . q{ -Xmx1024m -jar }
                          . $self->change_bam_header_jar()
                          . q{ VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true CREATE_MD5_FILE=true }
                          . qq{I=$old_bam O=$new_bam};
    if($sample_name){
       $change_header_cmd .= q{ SM="} . $self->_check_tag_value($sample_name) . q{"};
    }
    if($library_id){
       $change_header_cmd .= q{ LB="} . $self->_check_tag_value($library_id) . q{"};
    }
    if($study_name){
        $change_header_cmd .= q{ DS="} . $self->_check_tag_value($study_name);
        if($study_description){
           $change_header_cmd .= q{: }. $self->_check_tag_value($study_description);
        }
        $change_header_cmd .= q{"};
    }

    $self->log($change_header_cmd);
    my $change_rs =  system $change_header_cmd;
    if($change_rs){
      croak "change bam header failed: $change_header_cmd";
    }
}


sub strip_bam_tag {
    my $self = shift;

    my $old_bam = $self->working_dir().q{/}.$self->old_bam_filename();
    my $new_bam = $self->working_dir().q{/}.$self->bam_filename();

    my $strip_bam_tag_cmd = $self->java_cmd . q{ -Xmx1024m -jar }
                          . $self->strip_bam_tag_jar()
                          . q{ VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true CREATE_MD5_FILE=true }
                          . qq{I=$old_bam O=$new_bam }
                          .  q{STRIP=OQ};

    $self->log($strip_bam_tag_cmd);
    my $strip_rs =  system $strip_bam_tag_cmd;
    if($strip_rs){
      croak "Strip bam tag failed: $strip_bam_tag_cmd";
    }
}

sub _check_tag_value {
  my ($self, $tag_value) = @_;

  $tag_value =~ s/\n/\ /gmxs;
  $tag_value =~ s/\t/\ /gmxs;

  return $tag_value;
}


sub _get_number_of_reads_qc {
   my ($self, $bam) = @_;

   my $num_total_reads;

   my ($name, $path, $suffix) = fileparse($bam, q{.bam});

   my $bam_flagstats_file = File::Spec->catfile ( $self->qc_dir(), $name.q{_bam_flagstats.json} );

   if(-e $bam_flagstats_file ){

      my $flag_stats = from_json(slurp $bam_flagstats_file);
      $num_total_reads = $flag_stats->{num_total_reads};
   }

   return $num_total_reads;
}

sub convert_alignment_filter_json {
    my $self = shift;

    my $qc_in = $self->working_dir();
    my $qc_out   = $self->qc_dir();
    my $command = sprintf 'qc --check alignment_filter_metrics --id_run %i --position %i --qc_in %s --qc_out %s',
                               $self->id_run,
                               $self->lane,
                               $qc_in,
                               $qc_out;
    if (defined $self->tag_index) {
        $command .= ' --tag_index ' . $self->tag_index;
    }

    $self->log($command);
    system $command;
    if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
       croak "Alignment filter autoqc check failed: $ERRNO\n$command";
    }

}


no Moose;

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 NAME

npg_common::bam_align_irods

=head1 VERSION

$Revision$

=head1 SYNOPSIS

use npg_common::bam_align_irods;

my $bam_re_aln = npg_common::bam_align_irods->new_with_options();

$bam_re_aln->run();

=head1 DESCRIPTION

This moudule will download the irods file to the working directory, do the re-alignment and reload them back to irods

=head1 SUBROUTINES/METHODS

=head2 run

Main method. Control calling of other methods.

=head2 get_old_bam_from_irods

=head2 collate_old_bam

=head2 bam_realign

=head2 reload_bam

=head2 remove_obsolete_replication

=head2 change_bam_header_jar

=head2 change_bam_header

=head2 strip_bam_tag_jar

=head2 strip_bam_tag

=head2 non_consent_split

reload new aligned bam file to irods with bai file as well

=head2 delete_qc_bam_flagstats

=head2 reload_bam_flagstats

=head2 convert_alignment_filter_json

=head1 DIAGNOSTICS

=head1 DEPENDENCIES

=head1 CONFIGURATION AND ENVIRONMENT

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi, E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2011 GRL, by Guoying Qi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
