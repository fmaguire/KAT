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
    class KmerCounters{
    public:
        const char* hash_path;
        uint64_t nb_distinct;
        uint64_t nb_kept;

        KmerCounters(const char* _hash_path) :
            hash_path(_hash_path)
        {
            nb_distinct = 0;
            nb_kept = 0;
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

        // Final data (created by merging thread results)
        KmerCounters                    *kmer_counters;

    public:
        KmerFilter(KmerFilterArgs* _args) : args(_args)
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
            input_hash = NULL;

            // Create the final kmer counters
            kmer_counters = new KmerCounters(args->jellyfish_hash.c_str());

            if (args->verbose)
                cerr << "Sequence filter tool setup successfully." << endl;
        }

        ~KmerFilter()
        {
            if (kmer_counters)
                delete kmer_counters;

            kmer_counters = NULL;

            if (jfh)
                delete jfh;

            jfh = NULL;
        }

        void do_it()
        {

            if (args->verbose)
            {
                cerr << "Loading hashes..." << endl;
            }

            std::ostream* out_stream = args->verbose ? &cerr : (std::ostream*)0;

            // Load the hashes
            input_hash = jfh->loadHash(true, out_stream);


            if (args->verbose)
                cerr << endl
                     << "Hash loaded successfully." << endl
                     << "Processing...";

            // Setup iterator for this thread's chunk of the hash
            hash_query_t::iterator hashIt = input_hash->iterator_all();
            size_t hash_size = input_hash->get_size();
            uint_t kmer_len = input_hash->get_mer_len();
            uint_t val_len = input_hash->get_val_len();

            // Setup output hash, storage and dumper based on details from the input hash
            hash_reader_t* output_storage = NULL; //new inv_hash_storage_t(hash_size, 2*kmer_len,
                                    //val_len, input_hash->get_max_reprobe(), jellyfish::quadratic_reprobes);

            hash_writer_t output_hash(input_hash->get_nb_mers(), kmer_len, val_len, output_storage);

            // Go through this thread's slice for hash
            while (hashIt.next())
            {
                string kmer = hashIt.get_dna_str();
                uint64_t kmer_count = hashIt.get_val();

                if(keepKmer(kmer, kmer_count))
                {
                    uint64_t kmer_binary = jellyfish::parse_dna::mer_string_to_binary(kmer.c_str(), kmer_len);

                    output_hash.append(kmer_binary, kmer_count);
                }
            }

            if (args->verbose)
                cerr << "done." << endl;

            // Output hash to file
            ofstream_default hash_out_stream(args->output.c_str(), cout);
            output_hash.dump(&hash_out_stream);
            hash_out_stream.close();


            uint64_t nb_distinct = kmer_counters->nb_distinct;
            uint64_t nb_kept = kmer_counters->nb_kept;
            uint64_t nb_filtered = nb_distinct - nb_kept;

            *out_stream << "# Distinct Kmers Processed : " << nb_distinct << endl
                        << "# Discarded                : " << nb_filtered << endl
                        << "# Kept                     : " << nb_kept << endl;

        }


    private:


        bool keepKmer(string& kmer_seq, uint64_t kmer_count)
        {
            // Calculate GC%
            uint32_t gc_count = kat::gcCount(kmer_seq);

            // Are we within the limits
            bool in_gc_limits = args->low_gc <= gc_count && gc_count <= args->high_gc;

            // Check to see if we want to keep this sequence based on its GC content
            if ((in_gc_limits && args->discard) || (!in_gc_limits && !args->discard))
                return false;

            // Are we within the limits
            bool in_cvg_limits = args->low_count <= kmer_count && kmer_count <= args->high_count;

            // Return if we want to keep this sequence or not based on the coverage (and, if we've got this far, GC)
            return (in_cvg_limits && !args->discard) || (!in_cvg_limits && args->discard);
        }

    };
}
