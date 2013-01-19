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

/// \author Chris Saunders
///
#ifndef __BLOCKER_OPTIONS_HH
#define __BLOCKER_OPTIONS_HH


#include "print_double.hh"

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>


namespace FILTERTYPE {
    enum index_t {
        SITE,
        INDEL,
        BOTH
    };

    inline
    const char*
    label(const index_t x) {
        static const char* label[] = {"site","indel","locus"};
        return label[x];
    }

    inline
    const char*
    Label(const index_t x) {
        static const char* Label[] = {"Site","Indel","Locus"};
        return Label[x];
    }
}


struct FilterInfo {

    FilterInfo(const char* init_argname,
               const FILTERTYPE::index_t init_type,
               const char* init_label,
               const char* init_tag,
               const char* init_thresh,
               const bool init_is_max_thresh,
               const bool init_is_sample_value=false,
               const bool init_is_filter_if_missing =false)
        : argname(init_argname)
        , filter_type(init_type)
        , label(init_label)
        , tag(init_tag)
        , thresh(init_thresh)
        , is_max_thresh(init_is_max_thresh)
        , is_sample_value(init_is_sample_value)
        , is_filter_if_missing(init_is_filter_if_missing)
    {}

    std::string
    GetArgDescription() const {
        std::string s(is_max_thresh ? "Maximum " : "Minimum ");
        s += FILTERTYPE::label(filter_type);
        s += " " + tag;
        return s;
    }

    std::string argname; // used to enter a command-line argument
    FILTERTYPE::index_t filter_type;
    std::string label;
    std::string tag;
    print_double thresh;
    bool is_max_thresh; // filter is either a max or min allowed value
    bool is_sample_value; // if not sample value, then find the key in INFO for all samples
    bool is_filter_if_missing; // apply filter even if tag does not exist
};

std::ostream& operator<<(std::ostream& os,const FilterInfo& fi);



// gVCF nonvariant block settings, currently do not allow these to be set
// but putting them in blocker_opt makes this straightforward
//
struct NonvariantBlockOptions {
    NonvariantBlockOptions()
        : BlockFracTol("0.3")
        , BlockAbsTol(3)
        , BlockavgLabel("BLOCKAVG_min30p3a")
    {}

    print_double BlockFracTol;
    int BlockAbsTol;
    std::string BlockavgLabel;
};



struct BlockerOptions {

    BlockerOptions();

    ~BlockerOptions();

    void
    finalize_filters();

    bool
    is_chrom_depth() const {
        return (! ChromDepth.empty());
    }

    bool
    is_block_stats() const {
        return (! block_stats_file.empty());
    }

    enum filter_mode_t {
        FILTER_FULL,
        FILTER_NONE
    };

    std::ostream& outfp;
    bool is_skip_header;
    std::string max_chrom_depth_filter_tag;
    print_double max_chrom_depth_filter_factor;
    print_double min_nonref_blockable;
    std::string indel_conflict_label;
    std::string site_conflict_label;
    typedef std::map<std::string,double> cdmap_t;
    cdmap_t ChromDepth;
    std::string min_gqx;

    std::auto_ptr<FilterInfo> GQX_filter; // has to be a special case for now...
    std::vector<FilterInfo> filters;

    NonvariantBlockOptions nvopt;

    std::string block_stats_file;
    bool is_skip_blocks;
};


#endif
