#########
# Author:        gq1
# Maintainer:    $Author$
# Created:       2010 01 15
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

package npg_common::irods::Loader;

use strict;
use warnings;
use Moose;
use Carp;
use English qw(-no_match_vars);
use File::Basename;
use File::Spec;
use IPC::Open3;
use Digest::MD5 qw(md5);
use Perl6::Slurp;
use Symbol 'gensym';

with 'MooseX::Getopt';
with 'npg_common::roles::log';

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

Readonly::Scalar our $EXIT_CODE_SHIFT => 8;
Readonly::Scalar our $DEFAULT_RESOURCE   => q{};
Readonly::Scalar our $DEFAULT_ZONE       => q{seq};

## no critic (Documentation::RequirePodAtEnd ProhibitBacktickOperators)

=head1 NAME

npg_common::irods::Loader

=head1 VERSION

$LastChangedRevision$

=head1 SYNOPSIS

    my $loader = npg_common::irods::Loader->new({
       file         => $file_to_load,
       resource     => $destination_resourse,
       collection   => $destination_dir,
       new_filename => $destination_filename
       meta_data    => $meta_hashref,
    });

=head1 DESCRIPTION

This module load a file to irods server

you must run iinit command first

=head1 SUBROUTINES/METHODS

=head2 file

file to be uploaded, with full path

=cut

has 'file'        => (isa           => 'Str',
                      is            => 'rw',
                      required      => 0,
                      documentation => 'file to be uploaded, with path',
                     );
=head2 resource

which resource should the file be loaded to, otherwise use default

=cut

has 'resource'   =>  (isa           => 'Str',
                      is            => 'rw',
                      default       => $DEFAULT_RESOURCE,
                      documentation => 'resource to use',
                     );

=head2 zone

zone the file be loaded to, otherwise use default

=cut

has 'zone'       =>  (isa           => 'Str',
                      is            => 'rw',
                      default       => $DEFAULT_ZONE,
                      documentation => 'zone to use',
                     );

=head2 collection

which collection to load the file, otherwise current

=cut

has 'collection'  => (isa           => 'Str',
                      is            => 'rw',
                      required      => 0,
                      documentation => 'collection to load the file',
                     );

=head2 new_filename

given a new filename if you want to change it

=cut
has 'new_filename' =>(isa           => 'Str',
                      is            => 'rw',
                      required      => 0,
                      lazy_build    => 1,
                      documentation => 'new filename if you need change',
                     );

sub _build_new_filename{
   my $self = shift;
   return basename($self->file());
}

=head2 chmod_permissions

(array ref of) chmod permission strings, in the following format:
null|read|write|own userOrGroup, for example,
'null public' will give no permissions for public group or user

=cut

has 'chmod_permissions' => (is            => 'rw',
                            isa           => 'ArrayRef[Str]',
                            required      => 0,
                            default       => sub { [] },
                            documentation => 'chmod permission strings',
                         );

=head2 meta_data

a list of meta data as hash ref

=cut

has 'meta_data'   => (isa           => 'HashRef',
                      is            => 'rw',
                      required      => 0,
                      default       => sub { {} },
                      documentation => 'meta data to add',
                     );

=head2 _dry_run

dry_run flag

=cut

has '_dry_run'  => ( isa           => 'Bool',
                     is            => 'ro',
                     required      => 0,
                     documentation => 'dry_run flag',
                   );

=head2 run

main method to call, load the file and add meta data 

=cut

sub run{
  my $self = shift;

  my $des_file = $self->new_filename;

  my $des_dir = $self->collection();

  if($des_dir){

    $self->log('Making destination directory when neccesary');
    my $mkdir_cmd = qq{imkdir -p $des_dir};
    $self->_run_cmd($mkdir_cmd);

    $des_file = File::Spec->catfile($des_dir, $des_file);
  }

  $self->log('Checking md5 of the file');
  my $local_md5 = $self->_get_file_md5($self->file());
  my $remote_md5 = -d $self->file() ? undef : $self->_get_irods_md5($des_file);

  if( defined $remote_md5 && ( $remote_md5 eq $local_md5 ) ){
     $self->log( 'FILE EXISTS And MD5 CORRECT! NO NEED TO LOAD: ' . $self->file );
     $self->add_meta($des_file, $self->meta_data());
     return;
  }

  my $iput_cmd = q{iput};
  if (!-d $self->file()){ $iput_cmd .= q{ -K};} # -K is contraindicated for directories

  if($self->resource()){
     $iput_cmd .=q{ -R }.$self->resource();
  }

  $iput_cmd .= q{ -f};

  $iput_cmd .= q{ -r}; # required for directories

  $iput_cmd .= q{ } . $self->file();

  if( !-d $self->file() ){
     $iput_cmd .= qq{ $des_file};
  }elsif( $des_dir ){
     $iput_cmd .= qq{ $des_dir};
  }

  $self->_run_cmd($iput_cmd);

  if (! -d $self->file()) {
    $remote_md5 = $self->_get_irods_md5($des_file);
    if( !$remote_md5 || ($remote_md5 && $remote_md5 ne $local_md5) ) {
      croak "Loading failed, MD5 value not match.\nLocal file: $local_md5\niRODs file: $remote_md5";
    }else{
      $self->log('md5 checking correct after loading');
    }
  }

  if( my@permissions = @{$self->chmod_permissions()} ){
    eval {
       for my$permission (@permissions) {
         $self->run_set_permissions_command($permission, $des_file);
       }
      1;
    }  or do {
     carp "Problems with iRODS chmod: $EVAL_ERROR";
    };
  }

  $self->meta_data->{md5} = $local_md5;

  if(!-d $self->file()){
     $self->add_meta($des_file, $self->meta_data());
     $self->remove_obsolete_replication($des_file);
  }

  return 1;

}

=head2 add_meta

given a file on the server and a list of meta data, add the meta data to the file

=cut

sub add_meta{
  my ($self, $file, $meta_data) = @_;


  if( !$file ){

     my $collection = $self->collection();
     my $new_filename = $self->new_filename();
     $file = $collection ?  File::Spec->catfile($collection, $new_filename)
                         :  $new_filename;
  }

  if( !$meta_data ) {
     $meta_data = $self->meta_data();
  }

  my $remote_meta_data = $self->_check_meta_data($file);

  my %imeta_rmw_commands = ();
  my @imeta_add_commands = ();

  foreach my $att_name (keys %{$meta_data}){

     my $values = $meta_data->{$att_name};

     if( ref $values eq q{ARRAY} ){
       foreach my $value (@{$values}){
          my ($rmw_cmd, $add_cmd) = $self->_check_add_meta($file, $att_name, $value, $remote_meta_data);
          if($rmw_cmd){
            $imeta_rmw_commands{$rmw_cmd}++;
          }
          if($add_cmd){
            push @imeta_add_commands, $add_cmd;
          }
       }
     }else{
          my ($rmw_cmd, $add_cmd)  = $self->_check_add_meta($file, $att_name, $values, $remote_meta_data);
          if($rmw_cmd){
            $imeta_rmw_commands{$rmw_cmd}++;
          }
          if($add_cmd){
            push @imeta_add_commands, $add_cmd;
          }
     }

  }

  my $imeta_rmw_command_lines = join qq{\n}, (keys %imeta_rmw_commands, @imeta_add_commands);
  if (!$self->_dry_run) {
    $self->_run_imeta_commands($imeta_rmw_command_lines);
  } else {
    $self->log("The following imeta commands will be run:\n$imeta_rmw_command_lines");
  }

  return 1;
}

sub _run_imeta_commands{
  my ($self, $imeta_command_lines) = @_;

  my($wtr, $rdr, $err);
  $err = gensym;
  my $pid = open3($wtr, $rdr, $err, q{imeta} );

  print {$wtr} $imeta_command_lines.qq{\n} or croak q{Problems to pipe imeta sub command to imeta};

  close $wtr or croak q{Cannot close};

  waitpid $pid, 0;

  my $child_exit_status = $CHILD_ERROR >> $EXIT_CODE_SHIFT;

  $self->log(q{Error message:});
  while(my $line = <$err>){
    $self->log( $line );
  }

  $self->log(q{Output:});
  while(my $line = <$rdr>){
    $self->log( $line );
  }

  close $err or croak q{Cannot close};
  close $rdr or croak q{Cannot close};

  if($child_exit_status){
     croak 'Failed to run imeta command';
  }

  return 1;
}

sub _check_add_meta {
  my ($self, $file, $att_name, $value, $remote_meta_data) = @_;

  my $imeta_rmw_command_lines;
  my $imeta_add_command_lines;

  $value =~ s/^\s+//mxs;
  $value =~ s/\s+$//mxs;
  my $quoted = $value;
  if($quoted =~ /"/mxs && $quoted !~ /'/mxs){
    $quoted = qq{'$quoted'};
  }elsif($quoted =~ /'/mxs && $quoted !~ /"/mxs){
    $quoted = qq{"$quoted"};
  }else{
    $quoted =~ s/"/'/gmxs;
    $quoted = qq{"$quoted"};
  }

  my $add_meta_cmd = qq{add -d $file $att_name $quoted};

  if( exists $remote_meta_data->{$att_name} ){
     $imeta_rmw_command_lines = $self->_rm_metadata($file, $att_name);
     $imeta_add_command_lines =  $add_meta_cmd;
  }else{
     $imeta_add_command_lines =  $add_meta_cmd;
  }

  return ($imeta_rmw_command_lines, $imeta_add_command_lines);
}

=head2 get_permissions

Queries permissions for a file.
Returned a mapping of permissions to arrays of users and usergroups

=cut
sub get_permissions {
  my ($self, $filename) = @_;
  my $command = q{ils -A} .q{ }.$filename;
  my $pid = open3( undef, my $fh, undef, $command);
  my $error;
  my $permissions = {};
  if (<$fh>) {
    my $permissions_line = <$fh>;
    if (!$permissions_line) {
      $error = "Failed to get permissions from $command";
    } else {
      foreach (split /\s/smx, $permissions_line ) {
        next if !$_;
        my ($user, $perm) = $_ =~ /(\w+\#\w+):(\w+)/smx;
        if ($user && $perm) {
          push @{$permissions->{$perm}}, $user;
        }
      }
    }
  } else {
    $error =  "No input from $command";
  }
  waitpid $pid, 0;
  if( $CHILD_ERROR >> $EXIT_CODE_SHIFT){
     croak "Failed $command";
  }
  close $fh or croak "cannot close a handle to '$command' output: $ERRNO";
  if ($error) {
    croak $error;
  }
  return $permissions;
}

=head2 run_set_permissions_command

Sets permissions for a file

=cut
sub run_set_permissions_command {
  my ($self, $permissions, $filename) = @_;
  my $ichmod_cmd = q{ichmod }. $permissions.q{ }.$filename;
  $self->_run_cmd($ichmod_cmd);
  return;
}

=head2 restrict_file_permissions

remove all permissions from users and usergroups that are not the owners of the file

=cut
sub restrict_file_permissions {
  my ($self, $file) = @_;

  my $permissions = $self->get_permissions($file);
  my $users = {};
  foreach my $perm (keys %{$permissions}) {
    if ($perm ne 'own') {
      foreach my $u (@{$permissions->{$perm}}) {
        if (!exists $users->{$u}) {
          $users->{$u} = 1;
          $self->log("Will remove all permissions from $u for $file");
          if (!$self->_dry_run) {
            $self->run_set_permissions_command(qq{null $u}, $file);
          }
	}
      }
    }
  }
  return;
}

=head2 file_exists

Checks whether a file exists in the iRODS repository.

=cut
sub file_exists {
  my ($self, $file) = @_;
  if (!$file) {
    croak 'File to check should be given';
  }
  my $result = eval {
    $self->_run_cmd(q{ils } . $file);
    1;
  };
  return $result;
}

sub _run_cmd {
  my ($self, @cmd) = @_;
  $self->log(join q{ }, @cmd);
  my $output_meta = system @cmd;
  if($CHILD_ERROR){
     croak q{Failed: }.(join q{ }, @cmd).qq{\n$CHILD_ERROR $output_meta};
  }
  return;
}

sub _check_meta_data {
  my ($self, $des_file) = @_;

  my $imeta_ls_cmd = qq{imeta ls -d $des_file};

  my $pid = open3( undef, my $imeta_ls_out_fh, undef, $imeta_ls_cmd);

  my ($name, $value);
  my %meta_data = ();

  while (my $line = <$imeta_ls_out_fh> ) {

    chomp $line;
    #$self->log($line);
    if( $line =~ /^[-]{4}/mxs ){
       $name = undef;
       $value = undef;
    }elsif ( $line =~ /^attribute[:][ ](.+)/mxs ){
      $name = $1;
    }elsif( $line =~ /^value[:][ ](.+)/mxs ) {
      $value = $1;
      if(defined $name && defined $value) {
        $meta_data{$name}->{$value}++;
      }
    }
  }

  waitpid $pid, 0;

  if( $CHILD_ERROR >> $EXIT_CODE_SHIFT){
     croak "Failed: $imeta_ls_cmd";
  }

  close $imeta_ls_out_fh or croak "can not close imeta ls command output: $ERRNO";

  return \%meta_data;
}

sub _get_irods_md5 {
  my ($self, $des_file) = @_;

  #iRODs won't generate md5 if file size is zero
  #this loader will fail, because no md5 returned after iput command

  my %des_file_md5_hash = ();

  my $ils_cmd = qq{ils -L $des_file};

  my $pid = open3( undef, my $ils_out_fh, undef, $ils_cmd);

  while (my $line = <$ils_out_fh> ) {

    chomp $line;

    if ( $line =~ /ERROR[:][ ]lsUtil[:][ ]/mxs ){
       #file not exists in the server
       last;
    }elsif ( $line =~ /[ ]&[ ]/mxs ){
       #file exists in the server and md5 available
       $line =  <$ils_out_fh>;
       my ($md5) = $line =~ /^[ ]{4}(\w+)[ ]/mxs;
       $des_file_md5_hash{$md5}++;
    }
  }

  waitpid $pid, 0;

  if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     carp "$des_file not exist on irods";
  }

  close $ils_out_fh or croak "can not close ils command output: $ERRNO";

  my $num_md5_values = scalar keys %des_file_md5_hash;

  if( $num_md5_values > 1){
     #different md5 on different servers
     return;
  }elsif( $num_md5_values == 0 ){
     #file not exists or md5 not available
     return;
  }

  my @md5_values = keys %des_file_md5_hash;
  return shift @md5_values;
}

sub _get_file_md5 {
  my ($self, $file) = @_;

  my $md5_file = $file.q{.md5};
  if(-e $md5_file){
    my @md5_lines = slurp $md5_file, { chomp=>1 };
    my $md5_from_file = shift @md5_lines;
    if($md5_from_file){
       return $md5_from_file;
    }
  }

  open my $fh, '<', $file or croak "Can't open '$file': $ERRNO";
  binmode $fh;

  my $md5 = -d $file ? 0 : Digest::MD5->new->addfile($fh)->hexdigest;
  close $fh or croak "Can't close '$file': $ERRNO";

  return $md5;
}

=head2 get_collection_file_list

given a collection, return a hashref with a list of files in this collection

=cut

sub get_collection_file_list {
  my ( $self, $collection) = @_;

  my %file_list = ();

  my $ils_cmd = qq{ils $collection};
  $self->log($ils_cmd);

  my $pid = open3( undef, my $ils_out_fh, undef, $ils_cmd);

  while (my $line = <$ils_out_fh> ) {

    chomp $line;

    if ( $line =~ /ERROR[:][ ]lsUtil[:][ ]/mxs ){
       #collection not exists in the server
       last;
    }elsif ($line !~ /$collection/mxs){
       #chomp $line;
       $line =~ s/\s+//mxs;
       $file_list{$line}++;
    }
  }

  waitpid $pid, 0;

  if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak "Failed: $ils_cmd";
  }

  close $ils_out_fh or croak "can not close ils command output: $ERRNO";

  return \%file_list;
}

=head2 imeta_query

given a hash of irods meta key value list, return a list of objects

=cut

sub imeta_query {
   my ($self, $query_hashref) = @_;

   my $imeta_qu_cmd = q{imeta qu -z }. $self->zone().q{ -d };

   my @key_values;
   foreach my $field (keys %{$query_hashref} ){
      push @key_values, qq{$field = '}.$query_hashref->{$field}.q{'};
   }
   $imeta_qu_cmd .= join q{ and }, @key_values;

   $self->log($imeta_qu_cmd);

   my @file_list;

   my $pid = open3( undef, my $imeta_qu_out_fh, undef, $imeta_qu_cmd);

   while (my $line = <$imeta_qu_out_fh> ) {

     chomp $line;
     if ($line =~/dataObj\:/mxs){
       #chomp $line;
       my ($file) = $line =~ /dataObj[:][ ](.+)/mxs;
       push @file_list, $file;
     }
   }

   waitpid $pid, 0;

   if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
     croak "Failed: $imeta_qu_cmd";
   }

   close $imeta_qu_out_fh or croak "can not close imeta command output: $ERRNO";

   return \@file_list;
}

sub _rm_metadata {
   my ($self, $irods_file, $att_name ) = @_;

   my $rm_command = q{rmw -d } . $irods_file . q{ } . $att_name . q{ %};

   return $rm_command;
}

=head2 rm_file

delete a file of collection

=cut
sub rm_file
{
	my ($self, $file) = @_;
	$file ||= $self->file();
	my $rm_cmd = "irm -r $file";
	$self->_run_cmd($rm_cmd);
	return 1;
}

=head2 remove_obsolete_replication

remove any obsolete replication on irods server

=cut

sub remove_obsolete_replication {
   my ($self, $irods_file) = @_;

   my $rep_num_to_remove = $self->_get_rep_num_to_remove($irods_file);
   foreach my $rep_num (@{$rep_num_to_remove}){

   	 $self->remove_one_replication($irods_file, $rep_num);
   }
   return;
}

=head2 remove_one_replication

remove one replication on irods server

=cut

sub remove_one_replication {
   my ($self, $irods_file, $rep_num) = @_;

   my $irm_command = qq{irm -n $rep_num $irods_file};
   $self->log($irm_command);
   my $irm_rs =  system $irm_command;
   if($irm_rs){
          croak "Failed: $irm_command";
   }

   return 1;
}

sub _get_rep_num_to_remove {
    my ($self, $irods_file) = @_;

    my $replication_nums = $self->get_replication_numbers($irods_file);

    return $replication_nums->{obsolete_rep_nums};
}

=head2 get_replication_numbers

get replication number for a file in irods in two arrays, normal ones and obsolete ones

=cut
sub get_replication_numbers {
    my ($self, $irods_file) = @_;

    my @rep_nums = ();
    my @obsolete_rep_nums = ();

    my $ils_cmd = qq{ils -l $irods_file};

    my $pid = open3( undef, my $ils_out_fh, undef, $ils_cmd);

    while (my $line = <$ils_out_fh> ) {

       chomp $line;

       if ( $line =~ /ERROR[:][ ]lsUtil[:][ ]/mxs ){
          #file not exists in the server
          last;
       }else{

          my ($rep_num) = $line =~ /^[ ]{2}\w+\s+(\d+)/mxs;

          if(defined $rep_num){

              if( $line !~ /[ ]&[ ]/mxs ){
                push @obsolete_rep_nums, $rep_num;
              }else{
                push @rep_nums, $rep_num;
              }
          }
       }
    }

    waitpid $pid, 0;

    if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
       carp "Problems to list bam file $irods_file";
    }

    close $ils_out_fh or croak "can not close ils command output: $ERRNO";

    return {rep_nums => \@rep_nums, obsolete_rep_nums => \@obsolete_rep_nums};
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

=item  File::Basename

=item  File::Spec

=item  npg_common::roles::log

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
