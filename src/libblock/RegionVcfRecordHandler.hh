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

#ifndef __REGION_VCF_RECORD_HANDLER
#define __REGION_VCF_RECORD_HANDLER


#include "istream_line_splitter.hh"
#include "ref_util.hh"
#include "region_util.hh"
#include "VcfRecord.hh"

#include <iosfwd>
#include <string>



struct RegionVcfOptions {

    RegionVcfOptions();

    std::ostream& outfp;
    std::string ref_seq_file;
    region_util::region_t regions;
};



// base class for block-breaking tools:
//
struct RegionVcfRecordHandler {

    RegionVcfRecordHandler(const RegionVcfOptions& opt)
        : _opt(opt)
        , _scp(opt.ref_seq_file.c_str())
    {}

    virtual ~RegionVcfRecordHandler() {}

    void
    process_line(const istream_line_splitter& vparse);

private:

    virtual
    void
    process_block(const bool is_in_region,
                  const unsigned end,
                  VcfRecord& vcfr) const = 0;

    // does the current vcf record overlap with any regions?
    bool
    is_record_in_region(const istream_line_splitter& vparse);

    // if is_record_in_region is true, call this function repeatedly
    // to get the end position and the in/out region status of the
    // next intercepting intervals
    //
    // returns false when no more intervals exist
    //
    bool
    get_next_record_region_interval(bool& is_in_region,
                                    unsigned& end);

protected:
    const RegionVcfOptions& _opt;
    samtools_char_picker _scp;

private:
    std::string _last_chrom;
    bool _is_skip_chrom; // true when pos is past all regions in current chrom
    region_util::interval_group_t::const_iterator _rhead,_rend;

    unsigned _begin_pos,_end_pos; // used to provide the region intercept iterator
};


#endif
