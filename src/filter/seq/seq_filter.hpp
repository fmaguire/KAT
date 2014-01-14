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
#include <cstdio>
#include <vector>
#include <math.h>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include <jellyfish/hash.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include <seq_utils.hpp>

#include "seq_filter_args.hpp"

using std::fstream;
using std::ostringstream;
using std::cerr;
using std::ios;


namespace kat
{

    class SeqFilterCounter
    {
    public:
        uint64_t nb_records;
        uint64_t nb_kept;

        SeqFilterCounter() : nb_records(0), nb_kept(0)
        {
        }

        SeqFilterCounter(uint64_t _nb_records, uint64_t _nb_kept) :
            nb_records(_nb_records), nb_kept(_nb_kept)
        {
        }
    };

    class SeqFilter
    {
    private:
        // Arguments from user
        SeqFilterArgs *args;

        JellyfishHelper                 *jfh;
        hash_query_t                    *hash;		// Jellyfish hash


    public:
        SeqFilter(SeqFilterArgs* _args) : args(_args)
        {
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
            bool have_seq_1 = !args->seq_file_1.empty();
            bool have_seq_2 = !args->seq_file_2.empty();

            if (!have_seq_1 && !have_seq_2)
                throw "No sequence files were specified";

            if (!have_seq_1)
                throw "First sequence file not specified";

            std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

            // Load the jellyfish hash
            hash = jfh->loadHash(false, out_stream);

            // Filter the sequences and return the counts
            SeqFilterCounter counter = have_seq_2 ? processPairs() : processSingles();

            uint32_t nb_records_filtered = counter.nb_records - counter.nb_kept;

            *out_stream << "Processed : " << counter.nb_records << endl
                       << "Filtered  : " << nb_records_filtered << endl
                       << "Kept      : " << counter.nb_kept << endl;
        }


    private:

        SeqFilterCounter processSingles()
        {
            // Open file, create RecordReader and check all is well
            fstream in(args->seq_file_1.c_str(), ios::in);
            seqan::RecordReader<fstream, seqan::SinglePass<> > reader(in);

            // Create the AutoSeqStreamFormat object and guess the file format.
            seqan::AutoSeqStreamFormat formatTag;
            if (!guessStreamFormat(reader, formatTag))
            {
                stringstream ss;
                ss << "ERROR: Could not detect file format for: " << args->seq_file_1 << endl;
                throw ss.str();
            }

            seqan::CharString out_path = args->output_prefix.c_str();
            append(out_path, ".");
            append(out_path, getAutoSeqStreamFormatName(formatTag));


            /*seqan::SequenceStream hash_out(
                                    toCString(out_path),
                                    seqan::SequenceStream::WRITE,
                                    seqan::SequenceStream::FASTQ,
                                    seqan::SequenceStream::AUTO_TYPE);
            if (!isGood(hash_out))
            {
                stringstream ss;
                ss << "ERROR: Could not open the output file: " << toCString(out_path) << endl;
                throw ss.str();
            }*/

            FILE* hash_out = fopen(toCString(out_path), "w");

            uint64_t nb_records_read = 0;
            uint64_t nb_records_written = 0;

            std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

            *out_stream << "Processing sequences..." << endl;

            uint16_t k = hash->get_mer_len();

            // Processes sequences in batches of records to reduce memory requirements
            while(!atEnd(reader))
            {
                seqan::CharString id;
                seqan::Dna5String dna5seq;
                seqan::CharString qual;

                if (readRecord(id, dna5seq, qual, reader, formatTag) != 0)
                {
                    stringstream ss;
                    ss << "ERROR: Problem reading record from: " << args->seq_file_1 << endl;
                    throw ss.str();
                }

                nb_records_read++;

                string seq;
                convertSeq(dna5seq, seq);

                uint64_t seq_length = seq.length();
                uint64_t nb_kmers_seq = seq_length - k + 1;

                if (seq_length < k)
                {
                    cerr << seq << " from R1 is too short to compute coverage.  Sequence length is "
                         << seq_length << " and K-mer length is " << k << ". Discarding sequence." << endl;

                    continue;
                }

                // Check to see if we want to keep this sequence, and if so write to file
                bool keep_seq = keepSeq(seq);

                if (keep_seq)
                {
                    //if (writeRecord(hash_out, id, seq, qual) != 0)
                    //if (writeRecord(hash_out, id, seq, qual) != 0)
                    if (writeFastQRecord(hash_out, toCString(id), seq.c_str(), toCString(qual)) != 0)
                    {
                        stringstream ss;
                        ss << "ERROR: Problem writing record to: " << endl;
                        throw ss.str();
                    }

                    nb_records_written++;
                }
            }

            fclose(hash_out);

            return SeqFilterCounter(nb_records_read, nb_records_written);
        }

        SeqFilterCounter processPairs()
        {
            // Open file, create RecordReader and check all is well
            fstream in1(args->seq_file_1.c_str(), ios::in);
            seqan::RecordReader<fstream, seqan::SinglePass<> > reader1(in1);

            fstream in2(args->seq_file_2.c_str(), ios::in);
            seqan::RecordReader<fstream, seqan::SinglePass<> > reader2(in2);

            // Create the AutoSeqStreamFormat object and guess the file format.
            seqan::AutoSeqStreamFormat formatTag1;
            if (!guessStreamFormat(reader1, formatTag1))
            {
                stringstream ss;
                ss << "ERROR: Could not detect file format for: " << args->seq_file_1 << endl;
                throw ss.str();
            }

            // Create the AutoSeqStreamFormat object and guess the file format.
            seqan::AutoSeqStreamFormat formatTag2;
            if (!guessStreamFormat(reader2, formatTag2))
            {
                stringstream ss;
                ss << "ERROR: Could not detect file format for: " << args->seq_file_2 << endl;
                throw ss.str();
            }

            if (formatTag1.tagId != formatTag2.tagId)
            {
                stringstream ss;
                ss << "ERROR: Files are different format." << endl;
                throw ss.str();
            }

            // Make files names
            seqan::CharString out_path_1 = args->output_prefix.c_str();
            append(out_path_1, "_R1.");
            append(out_path_1, getAutoSeqStreamFormatName(formatTag1));

            seqan::CharString out_path_2 = args->output_prefix.c_str();
            append(out_path_2, "_R2.");
            append(out_path_2, getAutoSeqStreamFormatName(formatTag2));


            /*seqan::SequenceStream hash_out_1(
                                toCString(out_path_1),
                                seqan::SequenceStream::WRITE,
                                seqan::SequenceStream::FASTQ,
                                seqan::SequenceStream::AUTO_TYPE);
            if (!isGood(hash_out_1))
            {
                stringstream ss;
                ss << "ERROR: Could not open the output file: " << toCString(out_path_1) << endl;
                throw ss.str();
            }

            seqan::SequenceStream hash_out_2(
                                toCString(out_path_2),
                                seqan::SequenceStream::WRITE,
                                seqan::SequenceStream::FASTQ,
                                seqan::SequenceStream::AUTO_TYPE);
            if (!isGood(hash_out_2))
            {
                stringstream ss;
                ss << "ERROR: Could not open the output file: " << toCString(out_path_2) << endl;
                throw ss.str();
            }*/

            FILE* hash_out_1 = fopen(toCString(out_path_1), "w");
            FILE* hash_out_2 = fopen(toCString(out_path_2), "w");


            uint64_t nb_records_read = 0;
            uint64_t nb_records_written = 0;

            std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

            *out_stream << "Processing sequences..." << endl;

            // Processes sequences in batches of records to reduce memory requirements
            while(!atEnd(reader1) && !atEnd(reader2))
            {
                seqan::CharString id1;
                seqan::Dna5String dna5seq1;
                seqan::CharString qual1;

                seqan::CharString id2;
                seqan::Dna5String dna5seq2;
                seqan::CharString qual2;

                if (readRecord(id1, dna5seq1, qual1, reader1, formatTag1) != 0)
                {
                    stringstream ss;
                    ss << "ERROR: Problem reading record from: " << args->seq_file_1 << endl;
                    throw ss.str();
                }

                if (readRecord(id2, dna5seq2, qual2, reader2, formatTag2) != 0)
                {
                    stringstream ss;
                    ss << "ERROR: Problem reading record from: " << args->seq_file_2 << endl;
                    throw ss.str();
                }

                nb_records_read++;

                string seq1, seq2;
                convertSeq(dna5seq1, seq1);
                convertSeq(dna5seq2, seq2);


                uint16_t k = hash->get_mer_len();
                uint64_t seq1_length = seq1.length();
                uint64_t seq2_length = seq2.length();
                uint64_t nb_kmers_seq1 = seq1_length - k + 1;
                uint64_t nb_kmers_seq2 = seq2_length - k + 1;

                if (seq1_length < k)
                {
                    //cerr << seq1 << " from R1 is too short to compute coverage.  Sequence length is "
                    //     << seq1_length << " and K-mer length is " << k << ". Discarding sequence." << endl;

                    continue;
                }

                if (seq2_length < k)
                {
                    //cerr << seq2 << " from R2 is too short to compute coverage.  Sequence length is "
                    //    << seq2_length << " and K-mer length is " << k << ". Discarding sequence." << endl;

                    continue;
                }

                // Check to see if we want to keep this sequence, and if so write to file
                bool keep_seq = keepSeq(seq1) && keepSeq(seq2);

                if (keep_seq)
                {
                    if (writeFastQRecord(hash_out_1, toCString(id1), seq1.c_str(), toCString(qual1)) != 0)
                    //if (writeRecord(hash_out_1, id1, seq1, qual1) != 0)
                    //if (writeRecord(hash_out_1, id1, seq1) != 0)
                    {
                        stringstream ss;
                        ss << "ERROR: Problem writing record to: " << endl;
                        throw ss.str();
                    }

                    if (writeFastQRecord(hash_out_2, toCString(id2), seq2.c_str(), toCString(qual2)) != 0)
                    //if (writeRecord(hash_out_2, id2, seq2, qual2) != 0)
                    //if (writeRecord(hash_out_2, id2, seq2) != 0)
                    {
                        stringstream ss;
                        ss << "ERROR: Problem writing record from the provided sequence file." << endl;
                        throw ss.str();
                    }

                    nb_records_written++;
                }
            }

            fclose(hash_out_1);
            fclose(hash_out_2);

            return SeqFilterCounter(nb_records_read, nb_records_written);
        }

        // Hack because I can't figure out how to get SEQAN writeRecord to work properly for Fastq.
        int writeFastQRecord(FILE* out, const char* id, const char* seq, const char* qual)
        {
            fprintf(out, "@%s\n%s\n+\n%s\n", id, seq, qual);

            return 0;
        }


        // There's no substring functionality in SeqAn in this version (1.4.1).  So we'll just
        // use regular c++ string's for this bit.  The next version of SeqAn may offer substring
        // functionality, at which time I might change this code to make it run faster using
        // SeqAn's datastructures.
        void convertSeq(seqan::Dna5String& from, string& to)
        {
            stringstream ssSeq;
            ssSeq << from;
            to = ssSeq.str();
        }

        bool keepSeq(string& seq)
        {
            uint16_t k = hash->get_mer_len();
            uint64_t seq_length = seq.length();
            uint64_t nb_kmers = seq_length - k + 1;

            bool found = false;

            if (seq_length < k)
            {
                //cerr << seq << " is too short to compute coverage.  Sequence length is "
                //     << seq_length << " and K-mer length is " << k << ". Discarding sequence." << endl;

                return false;
            }
            else
            {
                for(uint64_t i = 0; i < nb_kmers; i++)
                {
                    string merstr = seq.substr(i, k);

                    // Jellyfish compacted hash does not support Ns so if we find one then just skip to the next mer
                    if (merstr.find("N") == string::npos)
                    {
                        const char* mer = merstr.c_str();
                        uint_t count = (*hash)[mer];

                        if (count > 0)
                        {
                            found = true;

                            if (args->discard)
                                return false;
                        }
                    }

                }
            }

            return (found && !args->discard) || (!found && args->discard);
        }

    };
}
