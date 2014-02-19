#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2010-11-15
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::run::BamExtraMeta;

use Moose;
use Carp;
use English qw(-no_match_vars);
use IPC::Open3;

use npg_common::irods::Loader;

with qw{
          MooseX::Getopt
          npg_common::roles::log
          npg_common::roles::software_location
       };

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $DEFAULT_ROOT_DIR   => q{/seq/};
Readonly::Scalar our $DEFAULT_LANE_NUMBERS   => 8;
Readonly::Scalar our $EXIT_CODE_SHIFT => 8;

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::run::BamExtraMeta

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head2 id_run

id_run
 
=cut
has 'id_run'       => (isa           => q{Int},
                       is            => q{rw},
                       documentation => q{id_run},
                      );

=head2 file_list

bam file list
 
=cut
has 'file_list'    => (isa           => q{ArrayRef},
                       is            => q{rw},
                       lazy_build    => 1,
                       documentation => q{bam file name list},
                      );

sub _build_file_list {
  my $self = shift;

  my @file_list = ();

  my $run_dir = $DEFAULT_ROOT_DIR.$self->id_run();

  my $ils_cmd = q{ils }.$run_dir;

  my $pid = open3( undef, my $ils_fh, undef, $ils_cmd );
  while (my $line = <$ils_fh> ){
    chomp $line;
    $line =~ s/\s+//mxs;
    if( $line =~ /bam$/mxs ){
       push @file_list, $run_dir.q{/}.$line;
    }
  }
   waitpid $pid, 0;
   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak qq{Failed: $ils_cmd};
   }

   close $ils_fh or croak "can not close ils output: $ERRNO";

  return \@file_list;
}

=head2 process

main method to call add meta data for one run in irods 
 
=cut

sub process{
  my $self = shift;

  my $id_run = $self->id_run();

  $self->log("Adding meta data for bam file for run $id_run");

  my $file_list = $self->file_list();
  $self->log('There are '.( scalar @{ $file_list} ). ' bam files to add');

  my $file_count = 0;
  foreach my $file (@{ $file_list}){

    $self->log("---- Processing file $file");

    my $meta = {};
    my $reference = $self->get_bam_reference($file);
    my $alignment = 0;
    if( $reference ) {
       $alignment = 1;
       $meta->{reference} = $reference;
       $self->log("Reference: $reference");
    }
    $meta->{alignment} = $alignment;
    my $loader = npg_common::irods::Loader->new({
       file        => $file,
       meta_data   => $meta,
       unique_meta_data => {reference => 1, alignment => 1},
    });
    $loader->add_meta($file, $meta);
    $file_count++;
  }

  $self->log("---- $file_count bam files processed for run $id_run");

  return 1;
}


sub _get_ref_from_bwa_pg {
   my ($self, $bwa_pg_line) = @_;

   my @ref = grep {/\/nfs\/\S+\/references/mxs} split /\s/mxs, $bwa_pg_line;
   if(scalar @ref){
      return $ref[0];
   }
   return;
}

sub _grep_bam_bwa_pg_line {
   my ( $self, $sam_header_fh ) = @_;

   while( my $line = <$sam_header_fh> ) {
      if( $line =~ /^\@PG/mxs && $line =~/ID\:bwa/mxs && $line !~ /ID\:bwa_sam/mxs) {
          chomp $line;
          return $line;
      }
   }
   return;
}

=head2 get_bam_reference

read bam header using samtools and check bwa or bwa_aln PG to get the reference used
 
=cut

sub get_bam_reference {
   my ( $self, $bam ) = @_;

   my $SAMTOOLS = $self->samtools_cmd();

   my $samtools_view_cmd = qq{iget $bam - | $SAMTOOLS view -H -};
   $self->log($samtools_view_cmd);

   my $pid = open3( undef, my $sam_header_fh, undef, $samtools_view_cmd);

   my $bwa_pg_line = $self->_grep_bam_bwa_pg_line($sam_header_fh);
   my $reference;
   if($bwa_pg_line){
      $reference = $self->_get_ref_from_bwa_pg( $bwa_pg_line );
   }

   waitpid $pid, 0;

   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak qq{Failed: $samtools_view_cmd};
   }

   close $sam_header_fh or croak "can not close samtools view output: $ERRNO";

   return $reference;
}

=head2 process_runlist

main method to call add meta data for one all runs in irods 
 
=cut

sub process_runlist {

  my $self = shift;
  my @run_list = ();

  my $ils_cmd = q{ils }.$DEFAULT_ROOT_DIR;

  my $pid = open3( undef, my $ils_fh, undef, $ils_cmd );
  while (my $line = <$ils_fh> ){
    chomp $line;
    my ($id_run) = $line =~ /\/seq\/(\d+)/mxs;
    if($id_run){
      push @run_list, $id_run;
      $self->log("Run: $id_run");
      npg_common::irods::run::BamExtraMeta->new(
                                           id_run           => $id_run,
                                           $self->resolved_paths,
                                          )->process();
    }
   }
   waitpid $pid, 0;
   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak qq{Failed: $ils_cmd};
   }
   close $ils_fh or croak "can not close ils output: $ERRNO";
   return \@run_list;
}

=head2 process_runlist_from_file

main method to call add meta data for a list of runs in irods 
 
=cut

sub process_runlist_from_file {
   my ($self, $list_file) = @_;

   open my $list_fh, '<', $list_file or croak "can not open the file $list_file";
   while(my $line = <$list_fh> ){
       chomp $line;
       npg_common::irods::run::BamExtraMeta->new (
                                            id_run  => $line,
                                            $self->resolved_paths,
                                                 )->process();
   }
   close $list_fh or croak "can not close file $list_file: $ERRNO"; ;
   return;
}
no Moose;

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item MooseX::Getopt

=item Carp

=item English -no_match_vars

=item Readonly

=item npg_common::irods::Loader

=item npg_common::roles::log

=item npg_common::roles::software_location

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2010 GRL, by Guoying Qi

This file is part of NPG.

NPG is free software: you can redistribute it and/or modify
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
