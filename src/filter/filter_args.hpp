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
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include <common_args.hpp>

#include "filter_main.hpp"

using std::string;
using std::cerr;
using std::cout;
using std::endl;

namespace kat
{
    const string KAT_FILTER_KMER_ID        = "kmer";
    const string KAT_FILTER_SEQ_ID         = "seq";

    const uint16_t MIN_ARGS = 0;

    class FilterArgs : public BaseArgs
    {
    private:
        string  mode_arg;
        int     mode_argc;
        char**  mode_argv;

    protected:

        // ***********************************************
        // These methods override BaseArgs virtual methods

        const string usage() const               { return "Usage: kat filter <mode>"; }
        const string shortDescription() const    { return "Filtering tools"; }
        const string longDescription() const
        {
            return  "First argument should be the filter mode you wish to use:\n" \
                    "  - kmer:    Filters a jellyfish kmer hash based on GC and kmer count limits\n" \
                    "  - seq:     Filters sequences from a file based on presence of a kmer in the sequence\n";
        }

        const string optionsDescription() const    { return ""; }

        vector<option>* longOptions()
        {
            vector<option>* long_options = new vector<option>();

            return long_options;
        }

        string shortOptions()                   { return ""; }
        void setOption(int c, string& option_arg) {}
        void processRemainingArgs(const vector<string>& remaining_args) {}
        const string currentStatus() const       { return ""; }

    public:

        // Default constructor
        FilterArgs() : BaseArgs(MIN_ARGS)
        {}

        // Constructor that parses command line options
        FilterArgs(int argc, char* argv[]) : BaseArgs(MIN_ARGS)
        {
            customParse(argc, argv);
        }

        string getMode() {
            return mode_arg;
        }

        int getModeArgC() {
            return mode_argc;
        }

        char** getModeArgV() {
            return mode_argv;
        }

        bool validMode(string mode_str)
        {
            return mode_str.compare(KAT_FILTER_KMER_ID) == 0 ||
                   mode_str.compare(KAT_FILTER_SEQ_ID) == 0;
        }


        void customParse(int argc, char *argv[])
        {
            if (argc <= 1)
            {
                error("No plot mode specified");
            }
            else if (validMode(string(argv[1]))) {

                mode_arg = argv[1];
                mode_argc = argc - 1;
                mode_argv = argv + 1;
            }
            else {

                // Let BaseArgs have a go, but make sure we fail after
                parse(argc, argv);

                error("Invalid command line arguments passed to \"kat plot\"");
                exit(1);
            }
        }
    };
}
