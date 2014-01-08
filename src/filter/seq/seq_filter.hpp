//  ********************************************************************
//  This file is part of KAT - the K-mer Analysis Toolkit.
//
//  KAT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  KAT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with KAT.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#pragma once

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <jellyfish/hash.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include <seq_utils.hpp>

#include "seq_filter_args.hpp"

using seqan::SequenceStream;
using seqan::StringSet;
using seqan::String;
using seqan::CharString;
using seqan::Dna5String;

namespace kat
{
    template<typename hash_t>
    class SeqFilter : public thread_exec
    {
    private:
        // Arguments from user
        SeqFilterArgs *args;

        JellyfishHelper                 *jfh;
        hash_t                          *hash;		// Jellyfish hash

        // Set as we go
        float low_gc_perc;
        float high_gc_perc;


    public:
        SeqFilter(SeqFilterArgs* _args) : args(_args)
        {
            // Some validation first
            if(args->high_count < args->low_count)
                throw "High kmer count value must be >= to low kmer count value";

            if(args->high_gc < args->low_gc)
                throw "High GC count value must be >= to low GC count value";

            if (args->verbose)
                cerr << "Setting up read filtering tool..." << endl;

            // Setup handles to load hashes
            jfh = new JellyfishHelper(args->jellyfish_hash);

            // Ensure hash is null at this stage (we'll load them later)
            hash = NULL;

            if (args->verbose)
                cerr << "Filter tool setup successfully." << endl;
        }

        ~SeqFilter()
        {
            if (jfh)
                delete jfh;

            jfh = NULL;
        }

        void do_it()
        {
            // Open file, create RecordReader and check all is well
            std::fstream in(args->seq_file_1.c_str(), std::ios::in);
            seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader1(in);

            std::fstream in(args->input_sequences.c_str(), std::ios::in);
            seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader2(in);



            // Create the AutoSeqStreamFormat object and guess the file format.
            seqan::AutoSeqStreamFormat formatTag1;
            if (!guessStreamFormat(reader1, formatTag))
            {
                std::cerr << "ERROR: Could not detect file format for: " << args->input_sequences << endl;
                return;
            }


            // Setup output streams for files
            if (args->verbose)
                *out_stream << endl;

            // Sequence K-mer counts output stream
            std::ostringstream count_path;
            count_path << args->output_prefix << "_counts.cvg";
            ofstream_default count_path_stream(count_path.str().c_str(), cout);

            // Average sequence coverage and GC% scores output stream
            std::ostringstream cvg_gc_path;
            cvg_gc_path << args->output_prefix << "_stats.csv";
            ofstream_default cvg_gc_stream(cvg_gc_path.str().c_str(), cout);
            cvg_gc_stream << "seq_name coverage gc% seq_length" << endl;

            CharString id1;
            Dna5String dna5seq1;
            CharString qual1;

            CharString id2;
            Dna5String dna5seq2;
            CharString qual2;

            uint32_t nb_records_read = 0;
            uint32_t nb_records_written = 0;

            if (args->verbose)
                *out_stream << "Processing sequences..." << endl;

            bool haveSeq1 = !args->seq_file_1.empty();
            bool haveSeq2 = !args->seq_file_2.empty();

            // Processes sequences in batches of records to reduce memory requirements
            while(!atEnd(reader1))
            {
                if (readRecord(id, dna5seq, qual, reader1, formatTag1) != 0)
                {
                    cerr << "ERROR: Problem reading record from the provided sequence file." << endl;
                    return;
                }

                if (readRecord(id, dna5seq, qual, reader1, formatTag1) != 0)
                {
                    cerr << "ERROR: Problem reading record from the provided sequence file." << endl;
                    return;
                }

                nb_records_read++;

                // There's no substring functionality in SeqAn in this version (1.4.1).  So we'll just
                // use regular c++ string's for this bit.  The next version of SeqAn may offer substring
                // functionality, at which time I might change this code to make it run faster using
                // SeqAn's datastructures.
                stringstream ssSeq;
                ssSeq << dna5seq;
                string seq = ssSeq.str();

                // Check to see if we want to keep this sequence, and if so write to file
                if (keepSeq(seq))
                {

                    if (seqan::writeRecord(out, id, seq, qual, formatTag) != 0)
                    {
                        cerr << "ERROR: Problem writing record from the provided sequence file." << endl;
                        return;
                    }

                    nb_records_written++;
                }
            }


            uint32_t nb_records_filtered = nb_records_read - nb_records_written;

            *out_stream << "Processed : " << nb_records_read << endl
                        << "Filtered  : " << nb_records_filtered << endl
                        << "Kept      : " << nb_records_written << endl;

        }


    private:


        bool keepSeq(string& seq)
        {

            if (args->gc)
            {

            }
            else
            // Calculate GC%
            float gc_perc = kat::calcGCPercentage(seq);

            // Are we within the limits
            bool in_gc_limits = args->low_gc <= gc_perc && gc_perc <= args->high_gc;

            // Check to see if we want to keep this sequence based on its GC content
            if ((in_gc_limits && args->discard) || (!in_gc_limits && !args->discard))
            {
                return false;
            }

            // Now calculate the mean coverage
            float mean_cvg = kat::calcMeanCoverage(seq, hash);

            // Convert to logscale
            float log_cvg = args->cvg_logscale ? log10(mean_cvg) : mean_cvg;

            // Are we within the limits
            bool in_cvg_limits = args->low_count <= log_cvg && log_cvg <= args->high_count;

            // Return if we want to keep this sequence or not based on the coverage (and, if we've got this far, GC)
            return (in_cvg_limits && !args->discard) || (!in_cvg_limits && args->discard);
        }

    };
}
