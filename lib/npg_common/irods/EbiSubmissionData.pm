#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2012-05-22
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::EbiSubmissionData;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use npg_common::irods::Loader;
use Date::Calc qw(:all);
use POSIX qw(strftime);
use DBI;

with qw{MooseX::Getopt npg_common::roles::log};

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::irods::EbiSubmissionData

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

my $submission = EbiSubmissionData->new();

$submission->process();

=head1 DESCRIPTION

Check ebi submission tracking database to get a list of bam files which status is public

and add submisison run accession number, submission accession number, submission file md5 and public date to irods meta data as ebi_run_acc, ebi_sub_acc, ebi_sub_md5 and ebi_sub_date

=head1 SUBROUTINES/METHODS

=head2 verbose

verbose

=cut

has 'verbose'           => (isa           => q{Bool},
                            is            => q{rw},
                            default       => 0,
                            documentation => q{verbose},
                           );

=head2 not_dry_run

not_dry_run

=cut

has 'not_dry_run'  => ( isa           => 'Bool',
                        is            => 'ro',
                        required      => 0,
                        documentation => 'imeta sub commands will be run and irods meta data will be updated',
                      );

=head2 from_date

which date to start checking

=cut

has 'from_date'         => (isa           => q{Str},
                            is            => q{rw},
                            lazy_build    => 1,
                            documentation => q{which date to start checking},
                           );
sub _build_from_date {

     my @today = Today();

     my @one_week_ago = Add_Delta_Days(@today, -7);## no critic (ProhibitMagicNumbers)

     return sprintf q{%4d-%02d-%02d}, @one_week_ago;
}

=head2 stop_date

which date to stop checking

=cut

has 'stop_date'         => (isa           => q{Str},
                            is            => q{rw},
                            lazy_build    => 1,
                            documentation => q{which date to stop checking},
                           );
sub _build_stop_date {

     my @today = Today();

     return sprintf q{%4d-%02d-%02d}, @today;
}

=head2 _dbh

database connection to submission database subtrack

=cut

has '_dbh'              => (isa           => q{DBI::db},
                            is            => q{ro},
                            lazy_build    => 1,
                            documentation => q{database connection to submission database subtrack},
                           );
sub _build__dbh {

   return DBI->connect(q{dbi:mysql:host=shap;port=3303;dbname=subtrack},
			              q{webuser_ro},
			              q{},
			              {
				             RaiseError => 1,
				             AutoCommit => 0,
			              }
			             );
}

=head2 _query_sql

sql to query database

=cut

has '_query_sql'         => (isa           => q{Str},
                            is            => q{ro},
                            lazy_build    => 1,
                            documentation => q{sql to query database},
                           );
sub _build__query_sql {

   return q{select s.run, s.lane, s.mux, s.file_name, s.ebi_run_acc, s.ebi_sub_acc, f.md5, f.public_date
from submission s
join sub_status stat on (stat.id = s.id and stat.is_current = 'Y')
join file f on (s.file_name = f.file_name)
where stat.status = 'P'
and s.file_name like '%.bam'
and f.public_date >= ?
and f.public_date <= ?};
}

=head2 process

main method to call to process

=cut
sub process{
  my $self = shift;

  $self->log('Check ebi submission tracking database now');

  $self->_lookup_submission_db();

  $self->_dbh()->disconnect();

  $self->log('Updating irods ebi submission data finished');

  return 1;
}

=head2 _lookup_submission_db 

lookup submission database

=cut
sub _lookup_submission_db {
    my $self = shift;

    my $sql = $self->_query_sql();

    my $sth = $self->_dbh->prepare($sql);

    $sth->execute($self->from_date(), $self->stop_date() ) or croak $sth->errstr;

    my $count_checked = 0;
    my $count_updated = 0;

    while (my $row = $sth->fetchrow_hashref() ) {
       $count_checked++;
        eval{
           my $imeta_commands = $self->_process_one_file($row);
           if(scalar @{$imeta_commands}){
              $count_updated++;
              if(!$self->not_dry_run()){
                 print join "\n", @{$imeta_commands} or croak 'cannot print';
                 print "\n" or croak 'cannot print';
              }
           }
           1;
        } or do {
           $self->log($EVAL_ERROR);
        };
   }

   $sth->finish();

   $self->log("$count_checked records checked and $count_updated records need to be updated");

   return;
}

=head2 _process_one_file

process one database row

=cut
sub _process_one_file{
   my ($self, $db_row) = @_;

   my $id_run    = $db_row->{run};
   my $lane      = $db_row->{lane};
   my $tag_index = $db_row->{mux};
   my $file_name = $db_row->{file_name};

   my $irods_file_name = $self->_get_irods_file($id_run, $lane, $tag_index, $file_name);
   $self->verbose && $self->log($irods_file_name);

   my $ebi_run_acc   = $db_row->{ebi_run_acc};
   my $ebi_sub_acc   = $db_row->{ebi_sub_acc};
   my $ebi_sub_md5   = $db_row->{md5};
   my $ebi_sub_date  = $db_row->{public_date};

   my $imeta_commands = [];

   if( !$ebi_run_acc || !$ebi_sub_acc || !$ebi_sub_md5 || !$ebi_sub_date ){
      croak "Not all ebi submission data available for $file_name";
   }else{
      my $ebi_submisison_data = {
                                  ebi_run_acc  => $ebi_run_acc,
                                  ebi_sub_acc  => $ebi_sub_acc,
                                  ebi_sub_md5  => $ebi_sub_md5,
                                  ebi_sub_date => $ebi_sub_date,
                                };
      my $current_irods_meta_data = $self->_get_current_irods_meta($irods_file_name);
      $imeta_commands = $self->get_imeta_commands($irods_file_name, $current_irods_meta_data, $ebi_submisison_data);

      if(scalar @{$imeta_commands}){
           if($self->not_dry_run()){
               my $imeta_commands_string = join qq{\n}, @{$imeta_commands};
               my $irods_loader = npg_common::irods::Loader->new();
               $irods_loader->_run_imeta_commands( $imeta_commands_string );
               $self->verbose && $self->log( qq{irods meta data updated for $irods_file_name!\n} );
           }
      }else{
          $self->verbose && $self->log( qq{irods meta data no need to be updated for $irods_file_name!} );
      }
   }

   return $imeta_commands;
}

=head2 _get_irods_file

check data from database and form the full irods file name

=cut
sub _get_irods_file {
   my ($self, $id_run, $lane, $tag_index, $file_name) = @_;

   if($file_name !~ /[.]bam$/mxs){
      croak "File $file_name is not a bam file";
   }

   if( $file_name !~ /^${id_run}_${lane}/mxs ){
      croak "File name $file_name doesn't match id_run $id_run or lane $lane";
   }

   if( defined $tag_index && $file_name !~ /\#$tag_index/mxs ){
      croak "Tag index $tag_index is wrong in file name $file_name";
   }

   if(!defined $tag_index && $file_name =~ /\#/mxs){
      croak "Tag index not given but exists in file name $file_name";
   }

   return "/seq/$id_run/$file_name";
}

=head2 _get_current_irods_meta

get current irods meta data

=cut
sub _get_current_irods_meta {
   my ($self, $irods_file) = @_;

   my $irods_loader = npg_common::irods::Loader->new();

   return $irods_loader->_check_meta_data($irods_file);
}

=head2 get_imeta_commands

comapre the new data and the meta data already in database to generate imeta commands

=cut
sub get_imeta_commands{
   my ($self, $irods_file, $current_irods_meta, $ebi_submisison_data) = @_;

   my @imeta_commands = ();

   foreach my $meta_name (sort keys %{$ebi_submisison_data}){

       my $new_meta_value = $ebi_submisison_data->{$meta_name};

       my $old_meta = $current_irods_meta->{$meta_name};

       if($old_meta && $old_meta->{$new_meta_value}){
          next;
       }elsif( $old_meta && scalar keys %{$old_meta} ){

          push @imeta_commands, qq{rmw -d $irods_file $meta_name "%"};

          my $old_values_string = (strftime q{[%Y-%m-%dT%H:%M:%S] }, localtime). (join q{,}, keys %{$old_meta});
          my $meta_name_history = $meta_name . q{_history};

          push @imeta_commands, qq{add -d $irods_file $meta_name_history "$old_values_string"};
       }
       push @imeta_commands, qq{add -d $irods_file $meta_name "$new_meta_value"};
   }

   if(scalar @imeta_commands){
      $self->verbose() && $self->log(":\n" . join "\n", @imeta_commands);
   }

   return \@imeta_commands;
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

=item Date::Calc

=item POSIX

=item DBI

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2012 GRL, by Guoying Qi

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
