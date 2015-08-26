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

#include <string.h>

#include "kmer/kmer_filter_main.hpp"
#include "seq/seq_filter_main.hpp"

#include "filter_args.hpp"
#include "filter_main.hpp"

using std::string;

using kat::FilterArgs;

// Start point
int kat::filterStart(int argc, char *argv[])
{
    // Parse args
    FilterArgs args(argc, argv);

    // Shortcut to mode
    string mode = args.getMode();

    // Pass remaining args to relevant child tool
    if (mode.compare(KAT_FILT_KMER_ID) == 0)
    {
        kat::kmerFilterStart(args.getModeArgC(), args.getModeArgV());
    }
    else if (mode.compare(KAT_FILT_SEQ_ID) == 0)
    {
        kat::seqFilterStart(args.getModeArgC(), args.getModeArgV());
    }

    return 0;
}
