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

#include <jellyfish/hash.hpp>
#include <jellyfish/parse_dna.hpp>
#include <jellyfish/compacted_hash.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/storage.hpp>
#include <jellyfish/jellyfish_helper.hpp>

#include <seq_utils.hpp>

#include "kmer_filter_args.hpp"

using jellyfish::storage_t;

namespace kat
{

    class Counter
    {
    public:
        uint64_t distinct;
        uint64_t total;

        Counter() : distinct(0), total(0)
        {}

        void increment(uint64_t total_inc)
        {
            distinct++;
            total += total_inc;
        }
    };

    class KmerFilter
    {
    private:
        // Arguments from user
        KmerFilterArgs *args;

        JellyfishHelper                 *jfh;
        hash_query_t                    *input_hash;		// Jellyfish hash
        hash_t                          output_hash;		// Jellyfish hash

    public:
        KmerFilter(KmerFilterArgs* _args) : args(_args)
        {
            // Some validation first
            if(args->high_count < args->low_count)
                throw "High kmer count value must be >= to low kmer count value";

            if(args->high_gc < args->low_gc)
                throw "High GC count value must be >= to low GC count value";

            if (args->verbose)
                cerr << "Setting up kmer filtering tool..." << endl;

            // Setup handles to load hashes
            jfh = new JellyfishHelper(args->jellyfish_hash);

            // Ensure hash is null at this stage (we'll load them later)
            input_hash = NULL;

            if (args->verbose)
                cerr << "done." << endl;
        }

        ~KmerFilter()
        {
            if (jfh)
                delete jfh;

            jfh = NULL;
        }

        void do_it()
        {

            if (args->verbose)
            {
                cerr << "Loading hash..." << endl;
            }

            std::ostream* out_stream = &cerr;

            // Load the hashes
            input_hash = jfh->loadHash(true, out_stream);


            if (args->verbose)
                cerr << endl
                     << "Hash loaded successfully." << endl;

            // Setup iterator for this thread's chunk of the hash
            hash_query_t::iterator hashIt = input_hash->iterator_all();
            size_t hash_size = input_hash->get_size();
            uint_t kmer_len = input_hash->get_mer_len();
            uint_t val_len = input_hash->get_val_len();
            uint_t max_reprobe_index = 62; //input_hash->get_max_reprobe(); // Can't get it working this way! :s
            uint_t nb_mers = input_hash->get_nb_mers();

            if(args->verbose)
            {
                cerr << "Attempting to create output hash(es) with the following settings: " << endl
                     << " mer length        = " << kmer_len << endl
                     << " val length        = " << val_len << endl
                     << " hash size         = " << hash_size << endl
                     << " max reprobe index = " << max_reprobe_index << endl
                     << " nb mers           = " << nb_mers << endl << endl;
            }

            // Setup output hash, storage and dumper based on details from the input hash
            inv_hash_storage_t storage1(hash_size, 2*kmer_len, val_len, max_reprobe_index, jellyfish::quadratic_reprobes);
            jellyfish::compacted_hash::writer<inv_hash_storage_t> hash_writer1(nb_mers, 2*kmer_len, val_len, &storage1);

            // Setup output hash, storage and dumper based on details from the input hash
            inv_hash_storage_t storage2(hash_size, 2*kmer_len, val_len, max_reprobe_index, jellyfish::quadratic_reprobes);
            jellyfish::compacted_hash::writer<inv_hash_storage_t> hash_writer2(nb_mers, 2*kmer_len, val_len, &storage2);

            if(args->verbose)
                cerr << "Successfully created writer and empty filtered hash.  Populating filtered hash...";

            // Setup counters
            Counter all_count, in_count, out_count;

            // Go through this thread's slice for hash
            while (hashIt.next())
            {
                string kmer = hashIt.get_dna_str();
                uint64_t kmer_count = hashIt.get_val();

                // Update stats
                all_count.increment(kmer_count);

                bool in_bounds = inBounds(kmer, kmer_count);

                in_bounds ? in_count.increment(kmer_count) : out_count.increment(kmer_count);

                if (!args->separate)
                {
                    // Workout if we want to keep this kmer or not
                    if((in_bounds && !args->invert) || (!in_bounds && args->invert))
                    {
                        hash_writer1.append(hashIt.get_key(), kmer_count);
                    }
                }
                else
                {
                    // We are just separating the kmers so update whichever hash is appropriate
                    if (in_bounds)
                    {
                        hash_writer1.append(hashIt.get_key(), kmer_count);
                    }
                    else
                    {
                        hash_writer2.append(hashIt.get_key(), kmer_count);
                    }
                }
            }

            if (args->verbose)
            {
                cerr << "done." << endl
                     << "Writing filtered hash to: " << args->output << " ... ";
            }

            ostringstream out1_name_stream;
            out1_name_stream << args->output << (args->separate ? ".in.jf" : "") << kmer_len << "_0";
            string out_file_1 = out1_name_stream.str();

            ostringstream out2_name_stream;
            out2_name_stream << args->output << (args->separate ? ".out.jf" : "") << kmer_len << "_0";
            string out_file_2 = out2_name_stream.str();


            // Output hash to file
            ofstream_default hash_out_stream(out_file_1.c_str(), cout);
            hash_writer1.write_header(&hash_out_stream);
            hash_writer1.dump(&hash_out_stream);
            hash_out_stream.close();

            if (args->separate)
            {
                // Output hash to file
                ofstream_default hash_out_stream2(out_file_2.c_str(), cout);
                hash_writer2.write_header(&hash_out_stream2);
                hash_writer2.dump(&hash_out_stream2);
                hash_out_stream2.close();
            }


            if (args->verbose)
                cerr << "done." << endl;

            *out_stream << "Distinct kmers in input      : " << all_count.distinct << endl
                        << "Distinct kmers in bounds     : " << in_count.distinct << endl
                        << "Distinct kmers out of bounds : " << out_count.distinct << endl
                        << "Total kmers in input         : " << all_count.total << endl
                        << "Total kmers in bounds        : " << in_count.total << endl
                        << "Total kmers out of bounds    : " << out_count.total << endl;
        }


    private:


        bool inBounds(string& kmer_seq, uint64_t kmer_count)
        {
            // Calculate GC%
            uint32_t gc_count = kat::gcCount(kmer_seq);

            // Are we within the limits
            bool in_gc_limits = args->low_gc <= gc_count && gc_count <= args->high_gc;

            // Are we within the limits
            bool in_cvg_limits = args->low_count <= kmer_count && kmer_count <= args->high_count;

            // Return if we want to keep this sequence or not based on the coverage (and, if we've got this far, GC)
            return in_gc_limits && in_cvg_limits;
        }

    };
}
