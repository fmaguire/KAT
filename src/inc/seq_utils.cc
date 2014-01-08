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

#include <seq_utils.hpp>

using std::string;


uint32_t kat::gcCount(string& seq)
{
    uint32_t seq_length = seq.length();
    uint32_t g_or_c = 0;

    for(uint32_t i = 0; i < seq_length; i++)
    {
        char c = seq[i];

        if (c == 'G' || c == 'g' || c == 'C' || c == 'c')
            g_or_c++;
    }

    return g_or_c;
}



// Calc GC%
float kat::calcGCPercentage(string& seq)
{
    uint32_t seq_length = seq.length();

    uint32_t g_or_c = 0;
    uint32_t a_or_t = 0;

    for(uint32_t i = 0; i < seq_length; i++)
    {
        char c = seq[i];

        // Only register a,t,g,cs, ignore everything else
        if (c == 'G' || c == 'g' || c == 'C' || c == 'c')
            g_or_c++;
        else if (c == 'A' || c == 'a' || c == 'T' || c == 't')
            a_or_t++;
    }

    return ((float)(g_or_c)) / ((float)(g_or_c + a_or_t));
}

// Calc mean coverage
float kat::calcMeanCoverage(string& seq, hash_t& hash)
{
    uint64_t seqLength = seq.length();
    uint64_t nbCounts = seqLength - kmer + 1;

    if (seqLength < kmer)
    {
        cerr << names[index] << ": " << seq << " is too short to compute coverage.  Sequence length is "
             << seqLength << " and K-mer length is " << kmer << ". Setting sequence coverage to 0." << endl;
        return -1.0;
    }
    else
    {
        uint64_t sum = 0;

        for(uint64_t i = 0; i < nbCounts; i++)
        {
            string merstr = seq.substr(i, kmer);

            // Jellyfish compacted hash does not support Ns so if we find one set this mer count to 0
            sum += merstr.find("N") != string::npos ? 0 : (*hash)[merstr.c_str()];
        }

        // Assumes simple mean calculation for sequence coverage for now... plug in Bernardo's method later.
        return (float)sum / (float)nbCounts;
    }

}
