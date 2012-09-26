// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

/// \file
///
/// \author Chris Saunders
///


#include "region_util.hh"

#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>

namespace region_util {

std::ostream& log_os(std::cerr);


static
void
parse_bedfile_regions(const std::string& region_file,
                      region_t& regions) {

    if (region_file.empty()) return;

    std::ifstream region_is(region_file.c_str());
    if (! region_is){
        log_os << "ERROR: Failed to open region file '" << region_file << "'\n";
        exit(EXIT_FAILURE);
    }

    unsigned line_no(0);
    bool is_parse_fail(false);

    std::string bed_chrom;
    unsigned bed_begin(0),bed_end(0);

    while(! region_is.eof()){
        ++line_no;

        region_is >> bed_chrom;
        if(region_is.fail()) {
            if(! region_is.eof()) is_parse_fail=true;
            break;
        }

        if(bed_chrom != "track" && bed_chrom != "browser") {
        
            region_is >> bed_begin >> bed_end;
            if(region_is.fail() || (bed_end<bed_begin)) {
                is_parse_fail=true;
                break;
            }
            
            if(regions.find(bed_chrom) == regions.end()) {
                regions[bed_chrom] = interval_group_t();
            }
            regions[bed_chrom].push_back(std::make_pair(bed_begin,bed_end));
        }

        region_is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if(is_parse_fail) {
        log_os << "ERROR: unexpected format in region bed file line no: " << line_no << "\n";
        exit(EXIT_FAILURE);
    }
}



/// sort and merge overlapping/adjacent regions:
static
void
merge_regions(region_t& regions) {

    interval_group_t tmp_regions;

    region_t::iterator chromi(regions.begin()), chromi_end(regions.end());
    for(;chromi!=chromi_end;++chromi) {
        interval_group_t& cregions(chromi->second);
        std::sort(cregions.begin(),cregions.end());

        const unsigned nr(cregions.size());
        for(unsigned i(0);i<nr;++i) {
            interval_t& ci(cregions[i]);
            interval_t& hi(i==0 ? ci : tmp_regions.back());
            if((i==0) || (ci.first > hi.second)) {
                tmp_regions.push_back(ci);
            } else if(ci.second > hi.second) {
                hi = std::make_pair(hi.first,ci.second);
            }
        }
        cregions.swap(tmp_regions);
        tmp_regions.clear();
    }
}



void
dump_regions(const region_t& regions,
             std::ostream& os) {

    region_t::const_iterator chromi(regions.begin()), chromi_end(regions.end());
    for(;chromi!=chromi_end;++chromi) {
        const std::string& chrom(chromi->first);
        const interval_group_t& cregions(chromi->second);

        const unsigned nr(cregions.size());
        for(unsigned i(0);i<nr;++i) {
            os << chrom << "\t" << cregions[i].first << "\t" << cregions[i].second << "\n";
        }
    }
}



void
get_regions(const std::string& region_file,
            region_t& regions) {

    parse_bedfile_regions(region_file,regions);
    merge_regions(regions);
}


}
