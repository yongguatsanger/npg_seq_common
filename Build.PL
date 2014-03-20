use strict;
use warnings;
use Module::Build;

my $class = Module::Build->subclass(code => <<'EOF');

  sub git_tag {
    my $version;
    my $gitver = q[./bin/gitver];
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
EOF

my  $requires = {
                'autodie'                         => '2.10',
                'base'                            => '2.12',
                'Carp'                            => '1.04',
                'Class::Std'                      => '0.0.8',
                'Cwd'                             => 0,
                'DBI'                             => '1.608',
                'Digest::MD5'                     => 0,
                'English'                         => '1.02',
                'Exporter'                        => '5.63',
                'Fatal'                           => 0,
                'Fcntl'                           => 0,
                'File::Basename'                  => 0,
                'File::Compare'                   => 0,
                'File::Copy'                      => 0,
                'File::Path'                      => 0,
                'File::Slurp'                     => '9999.13',
                'File::Spec'                      => 0,
                'File::Spec::Functions'           => 0,
                'File::stat'                      => 0,
                'File::Temp'                      => 0,
                'File::Which'                     => 0,
                'FindBin'                         => '1.47',
                'HTML::PullParser'                => '3.57',
                'IO::File'                        => 0,
                'IO::Tee'                         => 0,
                'IPC::Open3'                      => 0,
                'lib'                             => '0.5565',
                'List::Util'                      => '1.21',
                'Math::Round'                     => '0.06',
                'MIME::Lite'                      => '3.024',
                'MIME::Parser'                    => '5.427',
                'Moose'                           => '0.93',
                'Moose::Role'                     => '0.93',
                'Moose::Util'                     => 0,
                'Moose::Util::TypeConstraints'    => 0,
                'MooseX::ClassAttribute'          => 0,
                'MooseX::Getopt'                  => '0.39',
                'MooseX::Storage'                 => 0,
                'MooseX::StrictConstructor'       => '0.08',
                'Parallel::ForkManager'           => 0,
                'Perl6::Slurp'                    => '0.03',
                'POSIX'                           => '1.09',
                'Readonly'                        => '1.03',
                'Readonly::XS'                    => 0,
                'strict'                          => '1.03',
                'warnings'                        => '1.05',
                'YAML'                            => '0.68',

	              'npg::api::request'                       => 0,
                'npg::api::run'                           => 0,
          	    'npg::api::run_status_dict'               => 0,

	       	      'npg_qc::autoqc::results::bam_flagstats'  => 0,
                'npg_qc::Schema'                          => 0,

                'npg_tracking::data::reference'           => 0,
                'npg_tracking::data::reference::list'     => 0,
                'npg_tracking::data::reference::find'     => 0,
                'npg_tracking::data::reference::info'     => 0,
                'npg_tracking::glossary::lane'            => 0,
                'npg_tracking::glossary::run'             => 0,
                'npg_tracking::glossary::tag'             => 0,
                'npg_tracking::illumina::run::folder'     => 0,
                'npg_tracking::illumina::run::short_info' => 0,          
                'st::api::lims'                           => 0,
};

if ($^X !~ /^\/software\/perl-5\.8\.8/smx) {
  $requires->{'YAML::XS'} = 0;
}

my $builder =  Module::Build->new(

          'module_name'         => 'useful_modules',
          'dist_author'         => q(wtsi-npg <npg@sanger.ac.uk>),
          'dist_abstract'       => 'General purpose scripts and modules used by WTSI production pipeline and manual data processing',
          'license'             => 'gpl',
          'dist_version'        => $class->git_tag(),
          'script_files' => [glob q[scripts/*]],

          'configure_requires' => {
                'ExtUtils::CBuilder'              => 0,
                'Readonly'                        => 0,
                'Readonly::XS'                    => 0,
          },

          'build_requires' => {
                'ExtUtils::CBuilder'              => 0,
                'Test::Cmd'                       => 0,
                'Test::Compile'                   => 0,
                'Test::Distribution'              => 0,
                'Test::Deep'                      => '0.103',
                'Test::Exception'                 => '0.27',
                'Test::MockObject'                => 0,
                'Test::More'                      => '0.86',
                'Test::Pod'                       => 0,
                'Test::Pod::Coverage'             => 0,
                'Test::Perl::Critic'              => 0,
                'Test::Trap'                      => '0.2.0',
                'Test::Warn'                      => '0.11',
          },

          'requires'     => $requires,

          'dist'         => { COMPRESS => 'gzip', SUFFIX => 'gz', },
         );

$builder->create_build_script();

1;
