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

#include "BlockerVcfHeaderHandler.hh"

#include <cstring>

#include <iostream>



static
void
print_filter_header(const FilterInfo& filter,
                    const char* type,
                    std::ostream& os) {

    static const char* gl[] = { "less" , "greater" };

    os << "##FILTER=<ID=" << filter.label
       << ",Description=\"" << type << " " << filter.tag <<" is " << gl[filter.is_max_thresh]
       << " than " << filter.thresh.strval();
    if(filter.is_filter_if_missing) {
        os << " or not present";
    }
    os << "\">\n";
}



static
void
print_filter_header_set(const std::vector<FilterInfo>& filters,
                        std::ostream& os) {

    const unsigned nf(filters.size());
    for(unsigned i(0);i<nf;++i) {
        const FilterInfo& fi(filters[i]);
        const char* type(FILTERTYPE::Label(fi.filter_type));
        print_filter_header(fi,type,os);
    }
}



bool
BlockerVcfHeaderHandler::
is_skip_header_line(const istream_line_splitter& vparse) {

    // remove some header lines:
    const unsigned nw(vparse.n_word());
    const unsigned rms(_rmKeys.size());
    for(unsigned i(0);i<rms;++i) {
        for(unsigned word_index(0);word_index<nw;word_index++) {
            if(NULL != strstr(vparse.word[word_index],_rmKeys[i].c_str())) {
                return true;
            }
        }
    }
    return false;
}



void
BlockerVcfHeaderHandler::
process_final_header_line() {

    _os << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the region described in this record\">\n";
    _os << "##INFO=<ID=" << _opt.nvopt.BlockavgLabel 
        << ",Number=0,Type=Flag,Description=\"Non-variant site block."
        << " All sites in a block are constrained to be non-variant, have the same filter value,"
        << " and have all sample values in range [x,y] , y <= max(x+3,(x*(1+" << _opt.nvopt.BlockFracTol << ")))."
        << " All printed site block sample values are the minimum observed in the region spanned by the block\">\n";
            
    // new format tags:
    _os << "##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"RMS Mapping Quality\">\n";
    _os << "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">\n";

    // overlap tags:
    _os << "##FILTER=<ID=" << _opt.indel_conflict_label
           << ",Description=\"Locus is in region with conflicting indel calls.\">\n";
    _os << "##FILTER=<ID=" << _opt.site_conflict_label
        << ",Description=\"Site genotype conflicts with proximal indel call. This is typically a heterozygous SNV call made inside of a heterozygous deletion.\">\n";
            
    // special chrom-depth filter tag:
    if(_opt.is_chrom_depth()){
        _os << "##FILTER=<ID=" << _opt.max_chrom_depth_filter_tag
            << ",Description=\"Site depth is greater than " << _opt.max_chrom_depth_filter_factor.strval()
            << "x the mean chromosome depth\">\n";
        BlockerOptions::cdmap_t::const_iterator i(_opt.ChromDepth.begin()), i_end(_opt.ChromDepth.end());
        for(;i!=i_end;++i) {
            const std::string& chrom(i->first);
            const double chrom_thresh(i->second*_opt.max_chrom_depth_filter_factor.numval());
            _os << "##" <<  _opt.max_chrom_depth_filter_tag << "_" << chrom << "=" <<  chrom_thresh << "\n";
        }
    }
            
    // print the rest of the standard filter tags:
    if(NULL != _opt.GQX_filter.get()) {
        print_filter_header(*_opt.GQX_filter,"Locus",_os);
    }
    print_filter_header_set(_opt.filters,_os);
}
