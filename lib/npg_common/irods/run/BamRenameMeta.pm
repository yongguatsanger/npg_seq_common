#########
# Author:        gq1
# Maintainer:    $Author: gq1 $
# Created:       2010-11-15
# Last Modified: $Date: 2010-12-02 12:30:28 +0000 (Thu, 02 Dec 2010) $
# Id:            $Id: BamExtraMeta.pm 12048 2010-12-02 12:30:28Z gq1 $
# $HeadURL: svn+ssh://svn.internal.sanger.ac.uk/repos/svn/new-pipeline-dev/useful_modules/branches/prerelease-32.0/lib/npg_common/irods/run/BamExtraMeta.pm $
#

package npg_common::irods::run::BamRenameMeta;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use npg_common::irods::Loader;
use IPC::Open3;
use Perl6::Slurp;

with qw{MooseX::Getopt npg_common::roles::log};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision: 12048 $ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $DEFAULT_ROOT_DIR   => q{/seq/};
Readonly::Scalar our $EXIT_CODE_SHIFT => 8;

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::run::BamRenameMeta

=head1 VERSION

$LastChangedRevision: 12048 $

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head2 name

meta data name

=cut
has 'name'         => (isa           => q{Str},
                       is            => q{rw},
                       documentation => q{meta data name},
                      );

=head2 old_value

old meta value

=cut
has 'old_value'    => (isa           => q{Str},
                       is            => q{rw},
                       documentation => q{old meta data value},
                      );

=head2 new_value

new meta value

=cut
has 'new_value'    => (isa           => q{Str},
                       is            => q{rw},
                       documentation => q{new meta data value},
                      );

=head2 file_list

bam file list to rename meta data value

=cut
has 'file_list'    => (isa           => q{ArrayRef},
                       is            => q{rw},
                       lazy_build    => 1,
                       documentation => q{bam file name list},
                      );

sub _build_file_list {
  my $self = shift;

  my @file_list = ();

  my $meta_name = $self->name();
  my $old_value = $self->old_value();

  my $imeta_cmd = qq{imeta qu -z $npg_common::irods::Loader::DEFAULT_ZONE -d $meta_name = $old_value};
  $self->log($imeta_cmd);

  my $pid = open3( undef, my $imeta_fh, undef, $imeta_cmd );
  my $collection;
  my $obj;
  while (my $line = <$imeta_fh> ){
    chomp $line;
    if($line =~ /^collection\:/mxs){
       ($collection) = $line =~ /^collection\:\s*(.+)$/mxs;
    }elsif($line =~ /^dataObj\:/mxs){
       ($obj) = $line =~ /^dataObj\:\s*(.+)$/mxs;
       if($collection && $obj) {
          push @file_list, $collection.q{/}.$obj;
       }
    }
  }
  waitpid $pid, 0;
  if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak qq{Failed: $imeta_cmd};
  }

  close $imeta_fh or croak "can not close imeta output: $ERRNO";

  return \@file_list;
}

=head2 process

modify meta value

=cut

sub process{
  my $self = shift;

  my $name = $self->name();
  my $old_value = $self->old_value();
  my $new_value = $self->new_value();

  $self->log("Renaming $name meta data from $old_value to $new_value");

  my $file_list = $self->file_list();
  $self->log('There are '.( scalar @{ $file_list} ). ' bam files to rename');

  my $file_count = 0;
  foreach my $file (@{ $file_list}){

    $self->log("---- Processing file $file");
    # add new value for this meta data and keep the old one
    my $meta = {sample=>$new_value,};

    my $loader = npg_common::irods::Loader->new({
       file        => $file,
       meta_data   => $meta,
    });
    $loader->add_meta($file, $meta);
    $file_count++;
  }

  $self->log("---- $file_count bam files processed");

  return 1;
}

=head2 process_runlist

main method to call to read a list from csv file

=cut

sub process_runlist {
  my ($self, $file) = @_;
  my @file_contents = slurp $file, {irs=>"\n", chomp=>"\n"};
  foreach my $line (@file_contents){
     my @values = split m/,/msx, $line;
     my $old_value = shift @values;
     my $new_value = shift @values;
     if($old_value && $new_value){
        my $bam = npg_common::irods::run::BamRenameMeta->new(name=>'sample', old_value=>$old_value, new_value=>$new_value,);
        $bam->process();
     }
  }

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

=item npg_common::irods::Loader

=item npg_common::roles::log

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
