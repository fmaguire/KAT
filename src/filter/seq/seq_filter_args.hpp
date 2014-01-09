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
    const bool      DEFAULT_FILTER_SEQ_DISCARD       = false;
    const char*     DEFAULT_FILTER_SEQ_OUTPUT_PREFIX = "kat.filter.seq";
    const char*     DEFAULT_FILTER_SEQ_SEQ_FILE_1    = "";
    const char*     DEFAULT_FILTER_SEQ_SEQ_FILE_2    = "";


    const uint16_t  FILTER_SEQ_MIN_ARGS = 2;


    class SeqFilterArgs : public BaseArgs
    {
    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const
        {
            return "Usage: kat filter seq [options] <jellyfish_hash> <seq_file1> [<seq_file2>]";
        }

        const string shortDescription() const
        {
            return "Filters sequences based on whether a kmer is present in that sequence.";
        }

        const string longDescription() const
        {
            string long_desc = "The GCP tool allows the user to quickly and easily see content int a kmer hash that is " \
                    " differentiated by GC or kmer coverage.  Sometimes this allows the user to identify contamination, " \
                    " or otherwise suspcious / interesting content within a sample, and in that case the user may wish " \
                    " to either isolate or discard the kmers and/or external sequences associated with this content.  " \
                    " This tool \"ref\" allows the user to do exactly that.";

            return lineBreakString(long_desc, 78, "  ");
        }

        const string optionsDescription() const
        {
            ostringstream help_str;

            help_str << " -o, --output_prefix=path    Output prefix for the filtered file.  If paired files are input, the " \
                     << "                             suffix will start with \"_R?.\" (\"" << DEFAULT_FILTER_SEQ_OUTPUT_PREFIX << "\")" << endl
                     << " -d, --discard_selection     Discard sequences that contain kmers in the hash.  By default, this " \
                     << "                             tool will discard sequences that do not have kmers found in the hash." << endl;

            return help_str.str();
        }

        vector<option>* longOptions()
        {
            static struct option long_options_array[] =
            {
                {"output_prefix",       required_argument,  0, 'o'},
                {"discard_selection",   no_argument,        0, 'd'}
            };

            vector<option>* long_options = new vector<option>();

            for(uint8_t i = 0; i < 2; i++)
            {
                long_options->push_back(long_options_array[i]);
            }

            return long_options;
        }

        string shortOptions()
        {
            return "o:d:";
        }

        void setOption(int c, string& option_arg) {

            switch(c)
            {
            case 'o':
                output_prefix = string(optarg);
                break;
            case 'd':
                discard = true;
                break;
            }
        }

        void processRemainingArgs(const vector<string>& remaining_args)
        {
            jellyfish_hash = remaining_args[0];
            seq_file_1 = remaining_args[1];

            if (remaining_args.size() > 2)
                seq_file_2 = remaining_args[2];
        }

        const string currentStatus() const
        {
            ostringstream status;

            status  << "discard: " << discard << endl
                    << "output: " << output_prefix << endl
                    << "seq_file_1: " << seq_file_1 << endl
                    << "seq_file_2: " << seq_file_2 << endl
                    << "jellyfish_hash: " << jellyfish_hash << endl;

            return status.str().c_str();
        }

    public:

        bool            discard;
        string          output_prefix;
        string          seq_file_1;
        string          seq_file_2;
        string          jellyfish_hash;

        SeqFilterArgs() : BaseArgs(FILTER_SEQ_MIN_ARGS),
            discard(DEFAULT_FILTER_SEQ_DISCARD),
            output_prefix(DEFAULT_FILTER_SEQ_OUTPUT_PREFIX),
            seq_file_1(DEFAULT_FILTER_SEQ_SEQ_FILE_1),
            seq_file_2(DEFAULT_FILTER_SEQ_SEQ_FILE_2)
        { }

        SeqFilterArgs(int argc, char* argv[]) : BaseArgs(FILTER_SEQ_MIN_ARGS),
            discard(DEFAULT_FILTER_SEQ_DISCARD),
            output_prefix(DEFAULT_FILTER_SEQ_OUTPUT_PREFIX),
            seq_file_1(DEFAULT_FILTER_SEQ_SEQ_FILE_1),
            seq_file_2(DEFAULT_FILTER_SEQ_SEQ_FILE_2)
        { parse(argc, argv); }

        ~SeqFilterArgs() {}

    };
}




