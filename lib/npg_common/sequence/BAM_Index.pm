#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2010-10-05
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::sequence::BAM_Index;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use File::Basename;
use File::Spec::Functions qw(catfile);
use Cwd qw(abs_path);
use Readonly;

with qw /
         MooseX::Getopt
         npg_common::roles::log
         npg_common::roles::software_location
        /;

our $VERSION = '0';

Readonly::Scalar our $DEFAULT_COMMAND_OPTIONS => {
                               VALIDATION_STRINGENCY => 'SILENT',
                               VERBOSITY             => 'ERROR',
                                                 };

Readonly::Scalar our $BUILD_BAM_INDEX_JAR => qw[BuildBamIndex.jar];

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::sequence::BAM_Index

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

  my $bam = npg_common::sequence::BAM_Index->new({
                 input_bam     => 'input.bam',
               });
  $bam->process();
  $bam->output_bam();

=head1 DESCRIPTION

Picard BuildBamIndex to generate bam index file

=head1 SUBROUTINES/METHODS

=head2 java_Xmx

JAVA command Xmx value

=cut
has 'java_Xmx'       => ( isa             => 'Str',
                          is              => 'rw',
                          default         => '2000m',
                          documentation   => 'java command Xmx value',
                        );

=head2 jar_file

Picard BuildBamIndex jar file

=cut

has 'bam_index_jar_file' => (
            is            => 'rw',
            isa           => 'NpgCommonResolvedPathJarFile',
            default       => $BUILD_BAM_INDEX_JAR,
            coerce        => 1,
            documentation => 'Picard BuildBamIndex jar file',
);

=head2 bam_index_options

=cut

has 'bam_index_options'            => ( isa             => 'HashRef',
                                        is              => 'rw',
                                        default         => sub { $DEFAULT_COMMAND_OPTIONS },
                                        documentation   => 'Picard BuildBamIndex command options',
                                      );
=head2 input_bam

input bam file name with path
 
=cut

has 'input_bam'      => ( isa             => 'Str',
                          is              => 'rw',
                          required        => 1,
                          documentation   => 'input bam filename with path',
                        );

=head2 output_bam

output bam index file with path
 
=cut

has 'output_bam_index'  => ( isa             => 'Str',
                             is              => 'rw',
                             required        =>  0,
                             lazy_build      => 1,
                             documentation   => 'output bam index filename with path',
                            );
sub _build_output_bam_index {
   my $self = shift;

   my($filename, $directories, $suffix) = fileparse($self->input_bam(), qw(.bam));

   return $directories.$filename.q{.bai};
}

=head2 bam_index_cmd

construct the whole picard bam indexing command
 
=cut

sub bam_index_cmd {
   my $self = shift;

   my $cmd = $self->java_cmd.q{ -Xmx}.$self->java_Xmx();
   $cmd .= q{ -jar }.$self->bam_index_jar_file();
   $cmd .= q{ INPUT=}.$self->input_bam();

   my $options = $self->bam_index_options();
   foreach my $option_key (keys %{$options}){
       $cmd .= qq{ $option_key='}.$options->{$option_key}.q{'};
   }

   return $cmd;
}


=head2 process

main method to call

=cut

sub process {
  my $self = shift;

  if( ! -e $self->input_bam() ){
    croak 'Input bam file does not exist to build index: '.$self->input_bam();
  }

  my $bam_index_cmd = $self->bam_index_cmd();
  $self->log ( $bam_index_cmd );

  my $bam_index_rs = system $bam_index_cmd;
  if( $bam_index_rs != 0){
    croak "Picard Bam indexing failed:\n$bam_index_cmd";
  }

  if( ! -e $self->output_bam_index() ){
     croak 'output bam index file not generated: '.$self->output_bam_index();
  }

  $self->log('BAM index file generated: '.$self->output_bam_index());

  return 1;
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

=item File::Basename

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
