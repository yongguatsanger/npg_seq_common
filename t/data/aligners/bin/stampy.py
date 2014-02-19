#!/usr/bin/env bash

cat <<End-of-version-string
stampy v1.0.11 (r880), <gerton.lunter@well.ox.ac.uk>

Usage: ./stampy.py [options] [.fa files]


Option summary (--help for all):

Command options
 -G PREFIX file1.fa [...]           Build genome index PREFIX.stidx from fasta file(s) on command line
 -H PREFIX                          Build hash PREFIX.sthash
 -M FILE[,FILE]                     Map fastq file(s).  Use FILE.recaldata for recalibration if available
 -R FILE[,FILE]                     Compute recalibration data from fastq file(s)

Mapping/output options
 -g PREFIX                          Use genome index file PREFIX.stidx
 -h PREFIX                          Use hash file PREFIX.sthash
 -o FILE                            Write mapping output to FILE [stdout]
 --readgroup=ID:id,tag:value,...    Set read-group tags (ID,SM,LB,DS,PU,PI,CN,DT,PL)  (SAM format)
 --solexa, --solexaold, --sanger    Solexa read qualities (@-based); pre-v1.3 Solexa; and Sanger (!-based, default)
 --substitutionrate=F               Set substitution rate for mapping and simulation [0.001]
 --gapopen=N                        Gap open penalty (phred score) [40]
 --gapextend=N                      Gap extension penalty (phred score) [3]
 --bwaoptions=opts                  Options and <prefix> for BWA pre-mapper (quote multiple options)
 --bwamaxmismatch=N                 Max number of mismatches for BWA maps; -1=auto [-1]
 --bwatmpdir=S                      Set directory for BWA temporary files

General options
 --help                             Full help
 -v N                               Set verbosity level (0-3) [2]


End-of-version-string

