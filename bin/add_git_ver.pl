#!/usr/bin/env perl
#########
# Author:       kt6 
# Maintainer:    $Author$ 
# Created:       2010-10-05
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );
use File::Find;

our $VERSION = '0';

sub git_tag {
    my $version;
    my $gitver = q[./../scripts/gitver];
    if (!-e $gitver) {
      warn "$gitver script not found";
      $version = q[unknown];
    }
    if (!-x $gitver) {
      warn "$gitver script is not executable";
      $version = q[unknown];
    }
    if (!$version) {
      $version = `$gitver`;
      $version =~ s/\s$//smxg;
    }
    return $version;
}

sub wanted {

# $File::Find::dir is the current directory name,
# $_ is the current filename within that directory
# $File::Find::name is the complete pathname to the file.
  
    my $gitver = git_tag();
    if (-e && -f) {
      warn "Changing version of $_ to $gitver\n";
      my $backup = '.original';
      local $^I = $backup;
      local @ARGV = ($_);
      while (<>) {
        s/(\$VERSION\s*=\s*)('?\S+'?)\s*;/${1}'$gitver';/;
        print;
      }
      #unlink "$_$backup";
    } else {
      warn "File $_ not found\n";
    }
}

find(\&wanted, "./lib/test.pm");
exit 0;

__END__

=head1 NAME

add_git_ver.pl

=head1 VERSION

$LastChangedRevision$

=head1 USAGE

 add_git_ver 

=head1 CONFIGURATION

=head1 SYNOPSIS

=head1 DESCRIPTION

Add git version number from bin/git_ver to $VERSION and pod VERSION.

Source file could have

our $VERSION = '49.0-8-g742cd6d-dirty';

or 

use Readonly; Readonly::Scalar our $VERSION => do { my ($r) = q$Revision$ =~ /(\d+)/mxs; $r; };

POD has $LastChangedRevision$

=head1 VERSION

$LastChangedRevision$

=head1 SUBROUTINES/METHODS

=head1 REQUIRED ARGUMENTS

=head1 OPTIONS

=head1 EXIT STATUS

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item FindBin

=item File::Find

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Kate Taylor

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2014 GRL, by Kate Taylor

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
