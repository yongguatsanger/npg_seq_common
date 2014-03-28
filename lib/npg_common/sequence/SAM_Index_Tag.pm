#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2010-08-03
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::sequence::SAM_Index_Tag;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);

with 'MooseX::Getopt';

our $VERSION = '0';
## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::sequence::SAM_Index_Tag

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

=head1 DESCRIPTION

This modules add indexing sequence tag and quality score to first read in the sam file

=head1 SUBROUTINES/METHODS

=head2 sam

input sam file to add indexing tag sequence in

if missing, stdin will be used
 
=cut

has 'sam'             => (isa             => 'Str',
                          is              => 'rw',
                          default         => q{},
                          documentation   => 'Input sam to add indexing tag sequence in, if missing stdin will be used',
                         );

=head2 index_fastq_file

index fastq file name with path, the order of reads must be the same with the order in sam file

=cut

has 'index_fastq_file'         => (isa            => 'Str',
                                   is             => 'rw',
                                   documentation  => 'index fastq file name with path',
                                   required       =>  1,
                                  );

=head2 process

the main method to call

=cut

sub process{
   my $self = shift;

   #get index fastq file
   my $fastq_t = $self->index_fastq_file();
   my $fastq_t_fh;
   if( $fastq_t && (! -e $fastq_t ) ){
      croak "the given index fastq file $fastq_t not exists!";
   }elsif( $fastq_t ){
     open $fastq_t_fh, q{<}, $fastq_t or croak "Can not open file $fastq_t"; ## no critic (InputOutput::RequireBriefOpen)
   }


   #check where to get sam file and pre-process sam file
   my $sam_file = $self->sam();

   my $sam_fh;
   my $sam_line;

   #sam from a file or stdin
   if($sam_file){
     open $sam_fh, q{<}, $sam_file or croak "Can not open file:$sam_file"; ## no critic (InputOutput::RequireBriefOpen)
   }else{
     $sam_fh =q{STDIN};
   }

   #jump header section of sam
   $sam_line = <$sam_fh>;
   while($sam_line && $sam_line =~ /^@/mxs){
         print $sam_line or croak 'Can not print';
         $sam_line = <$sam_fh>;
   }

   #check other sam line
   while($sam_line){
      chomp $sam_line;
      print $self->_add_index_tag($sam_line, $fastq_t_fh ), "\n" or croak "Cannot print $sam_line";
      $sam_line = <$sam_fh>;
   }

   if( $sam_fh ){
     close $sam_fh or croak "Can not close file: $sam_file";
   }
   if( $fastq_t_fh ) {
     close $fastq_t_fh or croak "Can not close file: $fastq_t";
   }

   return;
}

=head2 _add_index_tag

add index tag sequence for one sam line

=cut

sub _add_index_tag {
   my ($self, $sam_line, $fastq_t_fh) = @_;

   chomp $sam_line;

   if(!$sam_line){
     croak 'The given sam line is empty';
   }

   my @sam_fields = split /\t/mxs, $sam_line;
   if( scalar  @sam_fields < 4 ) {## no critic (ValuesAndExpressions::ProhibitMagicNumbers)
     croak "SAM format wrong: $sam_line";
   }

   my $id_sam = shift @sam_fields;
   my $flag = shift @sam_fields;
   my $first_read_sam = ($flag & 64) >> 6; ## no critic (ValuesAndExpressions::ProhibitMagicNumbers)
   my $second_read_sam = ($flag & 128) >> 7; ## no critic (ValuesAndExpressions::ProhibitMagicNumbers)

   if($second_read_sam){

     return $sam_line;
   }

   my $indexing_id;
   my $indexing_seq;
   my $indexing_qul;

   if( $fastq_t_fh ){
     my $read_found = 0;
     while( !$read_found ){
         $indexing_id = <$fastq_t_fh>;
         if(!$indexing_id){
            carp 'index fastq file reach the end';
            last;
         }
         $indexing_seq = <$fastq_t_fh>;
         <$fastq_t_fh>;
         $indexing_qul = <$fastq_t_fh>;
         chomp $indexing_id;
         chomp $indexing_seq;
         chomp $indexing_qul;

        if($indexing_id !~ /^@/mxs || !$indexing_seq || !$indexing_qul){
           croak "Indexing fastq file format wrong: $indexing_id is not the first line of one read or no sequence or quality available";
         }

        if(length $indexing_seq != length $indexing_qul){
            croak "Indexing tag sequence length is not equal to quality length\n$indexing_seq\n$indexing_qul"
        }

        if($indexing_id =~ /^\@$id_sam/mxs){
          $read_found = 1;
        }
     }
     if(!$read_found){
        croak "Read id in sam file $id_sam does not have a indexing tag sequence in the indexing fastq file";
     }
   }

   if( $indexing_seq && $indexing_qul){
      $sam_line = $sam_line.qq{\tRT:Z:$indexing_seq\tQT:Z:$indexing_qul};
   }

   return $sam_line;
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
