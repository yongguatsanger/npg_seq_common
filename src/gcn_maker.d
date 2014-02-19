// Author:        John O'Brien jo3@sanger.ac.uk
// Maintainer:    $Author$
// Created:       2008-11-26
// Last Modified: $Date$
// Id:            $Id$
// $HeadURL$
// Original code by: Aylwyn Scally 2009


// Compile with:
// ~js10/gdc/bin/gdc -fversion=Posix gcn_maker.d -L ~js10/gdc/lib -lgtango -O3 -o gcn_maker

// Break a fasta reference file into windows, calculate the GC and N count for
// each window. Print to a file.


import tango.io.Console;
import tango.io.Stdout;
import tango.io.stream.TypedStream;
import tango.io.stream.LineStream;
import tango.io.stream.TextFileStream;
import tango.io.FileConduit;
import tango.util.ArgParser;
import Int = tango.text.convert.Integer;
import tango.text.Util;
import tango.text.Ascii;
import tango.math.Math;
import tango.util.collection.LinkSeq;
import tango.io.FilePath;
import tango.text.stream.LineIterator;


//bool pairfilt = false;
char[] na         = "NA";
char[] gcn_suffix = ".gcn";
char[] backup     = "backup_file";
int offset;

class RefGCN{
	TextFileInput input;
	char[] name;
	char[] windowsize;
	char[] gc;
	char[] n;
	int sequence_index;
	int sequence_position;

	this(char[] inputname){
		input = new TextFileInput(inputname);
	}

	bool getline(){
		char[] line;
		char[][] tok;
		if (input.readln(line) && line.length > 0){
			tok               = delimit(line, "\t");
			sequence_index    = Int.parse(tok[0]);
			name              = tok[1].dup;
			sequence_position = Int.parse(tok[2]);
			windowsize        = tok[3].dup;
			gc                = tok[4].dup;
			n                 = tok[5].dup;
			return true;
		}else
			return false;
	}
}

void main (char[][] args) {
	char[] gcn_outfile_name;
	char[] reference_file = null;
	int window_size       = 1000;
	int max_no_of_samples =    0;
	int posopt            =    1;
	bool window           = false;

	char[] usage = "usage: window_depth ref_file [-b=window] [-o=pos_offset] [-n=na_string] [-N=max_samples]";

	ArgParser parser = new ArgParser((char[] value, uint ordinal){
		reference_file = value;
	});
	parser.bind("-", "n=", (char[] value){na = value;});
	parser.bind("-", "b=", (char[] value){window_size = Int.parse(value); window = true;});
	parser.bind("-", "o=", (char[] value){posopt = Int.parse(value);});
	parser.bind("-", "N=", (char[] value){max_no_of_samples = Int.parse(value);});

	if (args.length < 2) {
		Stdout(usage).newline;
		return;
	}
	try
		parser.parse(args[1..$]);
	catch (Exception e){
		Stderr(usage).newline;
		return;
	}

	char[10] buf;
	auto windowtxt = Int.format(buf, window_size);
	gcn_outfile_name = reference_file ~ gcn_suffix ~ windowtxt;

	char[] name;
	int current_sample_base_count =  0;
	int in_sequence_marker        = -1;
	int sequence_position         =  0;
	int sample_gc_count           =  0;
	int sample_n_count            =  0;
	int sample_count              =  0;
	int sequence_index            = -1;
	int[char[]] reference_index;
	char[][] ls;

	auto path = new FilePath(gcn_outfile_name);

	TextFileOutput gcnout = new TextFileOutput(gcn_outfile_name);
	Stderr.formatln("Reading reference {}", reference_file);
	auto inp = new TypedInput!(char)(new FileConduit(reference_file));
	foreach (c; inp){

        // Found a new sequence in the reference fasta.
		if (c == '>'){

            // If we already have some data, print it to file.
			if (current_sample_base_count > 0){
				gcnout.formatln("{}\t{}\t{}\t{}\t{}\t{}", sequence_index, name, sequence_position, current_sample_base_count, sample_gc_count, sample_n_count);
				sample_count++;
				if (max_no_of_samples > 0 && sample_count >= max_no_of_samples)
					break;
			}

			sequence_index++;
			name.length               = 0;
			in_sequence_marker        = 0;
			sequence_position         = 0;
			sample_gc_count           = 0;
			sample_n_count            = 0;
			current_sample_base_count = 0;

        // If we're still reading a sequence header line, build the name.
		} else if (in_sequence_marker == 0){
            // At the end of the line, tidy up.
			if (c == '\n'){
				name = delimit(name, " ")[0]; // keep the first word
				reference_index[name.dup] = sequence_index;
				in_sequence_marker = 1;
			} else {
				name ~= c;
			}


		} else if (in_sequence_marker > 0 && !isSpace(c)){
			c = toUpper([c])[0];

			sample_gc_count += cast(int)(c == 'G' || c == 'C');
			sample_n_count  += cast(int)(c == 'N');
			sequence_position++;

			current_sample_base_count++;

            // When the sample reaches the window size, print it and start again.
			if (current_sample_base_count == window_size){
				gcnout.formatln("{}\t{}\t{}\t{}\t{}\t{}", sequence_index, name, sequence_position, current_sample_base_count, sample_gc_count, sample_n_count);

				current_sample_base_count = 0;
				sample_gc_count           = 0;
				sample_n_count            = 0;
				sample_count++;

				if (max_no_of_samples > 0 && sample_count >= max_no_of_samples)
					break;
			}
		}
	}

    // After parsing the reference fasta, print any leftover data.
	if (current_sample_base_count > 0){
		gcnout.formatln("{}\t{}\t{}\t{}\t{}\t{}", sequence_index, name, sequence_position, current_sample_base_count, sample_gc_count, sample_n_count);
	}
	gcnout.flush.close;
}
