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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <iostream>

#include <file_utils.hpp>

#include "seq_filter.hpp"
#include "seq_filter_args.hpp"
#include "seq_filter_main.hpp"

using std::cerr;

using kat::SeqFilterArgs;
using kat::SeqFilter;

// Start point
int kat::seqFilterStart(int argc, char *argv[])
{
    // Parse args
    SeqFilterArgs args(argc, argv);

    // Print command line args to stderr if requested
    if (args.verbose)
        args.print();

    // Make sure there is something to process
    if (!fileExists(args.jellyfish_hash))
    {
        cerr << endl << "Could not find jellyfish hash file at: " << args.jellyfish_hash << "; please check the path and try again." << endl << endl;
        return 1;
    }

    // Create the sequence coverage object
    SeqFilter<hash_query_t> ref(&args);

    // Do the work
    ref.do_it();

    // That's it!  Just exit now.
    return 0;

}
