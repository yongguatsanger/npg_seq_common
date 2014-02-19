#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2009-04-17
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::Hit_Sequence_SAM;

use strict;
use warnings;
use Carp;
use English qw(-no_match_vars);
use base qw(Class::Accessor);

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$LastChangedRevision$ =~ /(\d+)/mxs; $r; };

__PACKAGE__->mk_accessors(qw(filename hit_ids non_hit_ids num_sequences_total num_sequences_hit num_sequences_nonhit));

sub new {
  my ($class, $ref) = @_;
  $ref ||= {};
  bless $ref, $class;

  return $ref;
}

sub parsing_file{
  my ($self) = shift;

  my $hit_ids = $self->hit_ids();
  my $non_hit_ids = $self->non_hit_ids();

  if(!$hit_ids){

    $hit_ids = {};
    $self->hit_ids($hit_ids);
  }

  if(!$non_hit_ids){

    $non_hit_ids = {};
    $self->non_hit_ids($non_hit_ids);
  }

  my $num_sequences_hit = 0;
  my $num_sequences_nonhit = 0;

  my $filename = $self->filename();

  print {*STDERR} ">Reading sam file to get the aligned sequence ids\n" or croak $ERRNO;

  open my $sam_fh, q[<], $filename or croak "$filename:\n$ERRNO";## no critic (InputOutput::RequireBriefOpen) 

  while( my $line = <$sam_fh>) {

    chomp $line;

    #ignore header section
    if($line =~ /^@/mxs){
      next;
    }
    my @fields_line = split /\t/sxm, $line;

    my $sequence_id = $fields_line[0];
    my $rname =  $fields_line[2];

    if($rname eq q{*}){

      $num_sequences_nonhit++;

      if(!$hit_ids->{$sequence_id}){

        $non_hit_ids->{$sequence_id}++;
      }
    }else{

      $num_sequences_hit++;

      $hit_ids->{$sequence_id}++;
      delete $non_hit_ids->{$sequence_id};
    }
  }
  close $sam_fh or croak $ERRNO;

  my $num_sequences_total = $num_sequences_hit + $num_sequences_nonhit;
  $self->num_sequences_hit($num_sequences_hit);
  $self->num_sequences_nonhit($num_sequences_nonhit);
  $self->num_sequences_total($num_sequences_total);

  return 1;
}

1;

__END__

=head1 NAME

npg_common::Hit_Sequence_SAM

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

  my $oHitSequenceSAM = npg_common::Hit_Sequence_SAM->new();

=head1 DESCRIPTION

Class to parse SAM sequence alignment file to get a list of hit short read ids and non-hit short read ids

=head1 SUBROUTINES/METHODS

=head2 new - constructor to create a object

=head2 parsing_file - get the name of SAM file from the object, parsing the file and geting two lists of sequence ids (hit and non hit), store them in the object

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

Carp
English
Class::Accessor

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi, E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2009 GRL, by Guoying Qi

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.

=cut
