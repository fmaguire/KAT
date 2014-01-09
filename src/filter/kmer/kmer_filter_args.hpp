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

#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <stdint.h>
#include <vector>

#include <common_args.hpp>
#include <str_utils.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ostringstream;

namespace kat
{
    const uint64_t  DEFAULT_FILTER_KMER_LOW_COUNT    = 0;
    const uint64_t  DEFAULT_FILTER_KMER_HIGH_COUNT   = 10000;
    const uint16_t  DEFAULT_FILTER_KMER_LOW_GC       = 0;
    const uint16_t  DEFAULT_FILTER_KMER_HIGH_GC      = 31;
    const bool      DEFAULT_FILTER_KMER_DISCARD      = false;
    const char*     DEFAULT_FILTER_KMER_OUTPUT       = "kat.ref";

    const uint16_t  FILTER_KMER_MIN_ARGS = 1;


    class KmerFilterArgs : public BaseArgs
    {
    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const
        {
            return "Usage: kat filter kmer [options] <jellyfish_hash>";
        }

        const string shortDescription() const
        {
            return "Filters kmers in a jellyfish hash based on gc and kmer coverage limits.";
        }

        const string longDescription() const
        {
            string long_desc = "The GCP tool allows the user to quickly and easily see content int a kmer hash that is " \
                    " differentiated by GC or kmer coverage.  Sometimes this allows the user to identify contamination, " \
                    " or otherwise suspcious / interesting content within a sample, and in that case the user may wish " \
                    " to either isolate or discard the kmers associated with this content.  " \
                    " This kmer filtering tool allows the user to do exactly that.";

            return lineBreakString(long_desc, 78, "  ");
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << " -o, --output_prefix=path    Output prefix for the filtered file and hash. (\"" << DEFAULT_FILTER_KMER_OUTPUT << "\")" << endl
                     << " -l, --low_count=uint64      Low kmer count limit (" << DEFAULT_FILTER_KMER_LOW_COUNT << ")" << endl
                     << " -h, --high_count=uint64     High kmer count limit (" << DEFAULT_FILTER_KMER_HIGH_COUNT << ")" << endl
                     << " -x, --low_gc=uint16         Low GC count limit. (\"" << DEFAULT_FILTER_KMER_LOW_GC << "\")" << endl
                     << " -y, --high_gc=uint16        High GC count limit. (\"" << DEFAULT_FILTER_KMER_HIGH_GC << "\")" << endl
                     << " -d, --discard_selection     Discard kmers and sequences equal to or within the limits rather than outside the limits." << endl;

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"output",              required_argument,  0, 'o'},
                {"low_count",           required_argument,  0, 'l'},
                {"high_count",          required_argument,  0, 'h'},
                {"low_gc",              required_argument,  0, 'x'},
                {"high_gc",             required_argument,  0, 'y'},
                {"discard_selection",   no_argument,        0, 'd'}
            };

            vector<option>* long_options = new vector<option>();

            for(uint8_t i = 0; i < 6; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return "o:l:h:x:y:d";
        }

        void setOption(int c, string& option_arg) {

            switch(c)
            {
            case 'o':
                output = string(optarg);
                break;
            case 'l':
                low_count = strToInt64(optarg);
                break;
            case 'h':
                high_count = strToInt64(optarg);
                break;
            case 'x':
                low_gc = strToInt16(optarg);
                break;
            case 'y':
                high_gc = strToInt16(optarg);
                break;
            case 'd':
                discard = true;
                break;
            }
        }

        void processRemainingArgs(const vector<string>& remaining_args)
        {
            jellyfish_hash = remaining_args[0];
        }

        const string currentStatus() const
        {
            ostringstream status;

            status  << "low_count: " << low_count << endl
                    << "high_count: " << high_count << endl
                    << "low_gc: " << low_gc << endl
                    << "high_gc: " << high_gc << endl
                    << "discard: " << discard << endl
                    << "output: " << output << endl
                    << "jellyfish_hash: " << jellyfish_hash << endl;

            return status.str().c_str();
        }

    public:

        uint64_t        low_count;
        uint64_t        high_count;
        uint32_t        low_gc;
        uint32_t        high_gc;
        bool            discard;
        string          output;
        string          jellyfish_hash;

        KmerFilterArgs() : BaseArgs(FILTER_KMER_MIN_ARGS),
            low_count(DEFAULT_FILTER_KMER_LOW_COUNT),
            high_count(DEFAULT_FILTER_KMER_HIGH_COUNT),
            low_gc(DEFAULT_FILTER_KMER_LOW_GC),
            high_gc(DEFAULT_FILTER_KMER_HIGH_GC),
            discard(DEFAULT_FILTER_KMER_DISCARD),
            output(DEFAULT_FILTER_KMER_OUTPUT)
        { }

        KmerFilterArgs(int argc, char* argv[]) : BaseArgs(FILTER_KMER_MIN_ARGS),
            low_count(DEFAULT_FILTER_KMER_LOW_COUNT),
            high_count(DEFAULT_FILTER_KMER_HIGH_COUNT),
            low_gc(DEFAULT_FILTER_KMER_LOW_GC),
            high_gc(DEFAULT_FILTER_KMER_HIGH_GC),
            discard(DEFAULT_FILTER_KMER_DISCARD),
            output(DEFAULT_FILTER_KMER_OUTPUT)
        { parse(argc, argv); }

        ~KmerFilterArgs() {}

    };
}




