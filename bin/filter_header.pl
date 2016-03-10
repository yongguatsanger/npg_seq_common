use strict;
use warnings;
use autodie;
use Getopt::Long;
use Carp;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );
use st::api::lims;

our $VERSION = '0';


my ($rpt,$run,$pos,$tag,$af1,$af2,$lims,$sample,$library,$study,);

&GetOptions("rpt:s"  => \$rpt,);

# Parse file name
if($rpt =~ /^(\d+)_(\d)(_(\w+))?(\#(\d+))?(_(\w+))?$/mxs){
    ($run,$pos,$tag,$af1,$af2) = ($1,$2,$6,$4,$8);
    if(defined $run && defined $pos && defined $tag){
        if($tag == 0){
            $lims = st::api::lims->new(id_run => $run, position => $pos);
            ($sample, $library, $study) = _get_limsm($lims);
        }else{
            $lims = st::api::lims->new(id_run => $run, position => $pos, tag_index => $tag);
            ($sample, $library, $study) = _get_limsi($lims);
        }
    }elsif(defined $run && defined $pos){
        $lims = st::api::lims->new(id_run => $run, position => $pos );
        ($sample, $library, $study) = _get_limsi($lims);
    }else{ die; }
}else{ die; }

# Process header stream
while(<>){
    chomp;
    if(/^\@RG/){
        my @l = split(/\t/);
        my $i;
        my ($sm, $lb, $ds);
        for($i = 0; $i < @l; $i++){
            if($l[$i] =~ /^SM:(.*)$/xms){
                $sm = $1;
                $l[$i] = q[SM:] . $sample;
                _compare_tag_values(q[SM],$sm, $sample, $rpt)
            }elsif($l[$i] =~ /^LB:(.*)$/xms){
                $lb = $1;
                $l[$i] = q[LB:] . $library;
                _compare_tag_values(q[LB], $lb, $library, $rpt)
            }elsif($l[$i] =~ /^DS:(.*)$/xms){
                $ds = $1;
                $l[$i] = q[DS:] . $study;
                _compare_tag_values(q[DS], $ds, $study, $rpt)
            }
        }
        print join(qq[\t], @l) ."\n";
    }else{
        print "$_\n";
    }
}

sub _get_limsm {
    my ($lims) = @_;
    my(@samples,@studies,%s,);
    foreach my $plex ($lims->children) {
        next if $plex->is_phix_spike;
        my ($sample_name,$library_id,$study) = _get_limsi($plex);
        push @samples, $sample_name;
        push @studies, $study if ! defined $s{$study};
        $s{$study}++;
    }
    my $sample_list = join ',', @samples;
    my $study_list  = join ',', @studies;
    return($sample_list, 'unknown', q[Study ]. $study_list);
}

sub _get_limsi {
    my ($lims) = @_;
    my $sample_name       = _check_tag_value($lims->sample_publishable_name());
    my $library_id        = _check_tag_value($lims->library_id());
    my $study_name        = _check_tag_value($lims->study_publishable_name());
    my $study_description = _check_tag_value($lims->study_description());
    if($lims->is_phix_spike){
        $study_description = 'SPIKED_CONTROL'
    }
    return($sample_name, $library_id, $study_name. q[: ].$study_description);
}

sub _check_tag_value {
    my ($tag_value) = @_;
    $tag_value =~ s/\n/\ /gmxs;
    $tag_value =~ s/\t/\ /gmxs;
    return $tag_value;
}

sub _compare_tag_values {
    my ($tag, $file_val, $lims_val, $rpt) = @_;
    if($file_val ne $lims_val){
        my $msg = qq[WARNING: LIMS data for: '$rpt', has changed: ];
        $msg .= qq[$tag old: '$file_val', ];
        $msg .= qq[$tag new: '$lims_val'.];
        carp q[$msg];
    }
}

__END__

=head1 NAME

filter_header.pl

=head1 VERSION

=head1 USAGE

samtools view -H 6551_1#1.cram | filter_header.pl -rpt 6551_1#1

=head1 CONFIGURATION

=head1 SYNOPSIS

=head1 DESCRIPTION

Reads in a BAM/CRAM file header and prints out an updated version of
the SM, LB and SM tags in the @RG section. Inteded to be used for
re-headering of such files.

Uses the value of the -rpt argument to obtain meta data from the
LIMS to update the values of the SM, LB and DS tags in the @RG
section. At the same time, it compares what is already there (or not)
and reports it by printing the results in STDERR. For the value of
-rpt to be valid it has to be of the form:

run_position[#tag[_split]|_split].

The input and output is the same: a BAM/CRAM file header in human-
readable format (e.g. using samtools view -H).

=head1 SUBROUTINES/METHODS

=over

=item  _get_limsm

Return multiple LIMS values as a concatenated list.

=item _get_limsi

Return individual LIMS values.

=item _check_tag_value

Remove '\t' and '\n' characters contained in LIMS information.

=item _compare_tag_values

Compare the value obtained from the LIMS vs the value of SM, LB and
DS present in the header and prints a message if they are different.

=back

=head1 REQUIRED ARGUMENTS

=head1 OPTIONS

=head1 EXIT STATUS

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item autodie

=item FindBin

=item lib

=item Getopt::Long

=item Carp

=item st::api::lims

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Carol Scott E<lt>ces@sanger.ac.ukE<gt> and Ruben Bautista E<lt>rb11@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015 GRL

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
