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

#include "BlockerOptions.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os,const FilterInfo& fi){

    os << "FilterInfo\n"
       << "\tlabel " << fi.label << "\n"
       << "\ttag " << fi.tag << "\n"
       << "\tthresh " << fi.thresh << "\n"
       << "\tis_max_thresh " << fi.is_max_thresh << "\n"
       << "\tis_sample_value " << fi.is_sample_value << "\n"
       << "\tis_filter_if_missing " << fi.is_filter_if_missing << "\n";
    return os;
}



BlockerOptions::
BlockerOptions()
    : outfp(std::cout)
    , is_skip_header(false)
    , max_chrom_depth_filter_tag("MaxDepth")
    , max_chrom_depth_filter_factor("3.0")
    , min_nonref_blockable("0.2")
    , indel_conflict_label("IndelConflict")
    , site_conflict_label("SiteConflict")
    , min_gqx("20.0")
    , is_skip_blocks(false)
{
    // set default filters:
    filters.push_back(FilterInfo("min-mq",FILTERTYPE::SITE,"LowMQ","MQ","20.0",false));
    filters.push_back(FilterInfo("min-qd",FILTERTYPE::BOTH,"LowQD","QD","3.73",false));
    filters.push_back(FilterInfo("max-site-fs",FILTERTYPE::SITE,"HighFS","FS","60.0",true));
    filters.push_back(FilterInfo("max-hapscore",FILTERTYPE::SITE,"HighHaplotypeScore","HaplotypeScore","13.0",true));
    filters.push_back(FilterInfo("min-mqrs",FILTERTYPE::SITE,"LowMQRankSum","MQRankSum","-12.5",false));
    filters.push_back(FilterInfo("min-site-rprs",FILTERTYPE::SITE,"LowReadPosRankSum","ReadPosRankSum","-2.386",false));
    
    filters.push_back(FilterInfo("max-indel-fs",FILTERTYPE::INDEL,"HighIndelFS","FS","200.0",true));
    filters.push_back(FilterInfo("min-indel-rprs",FILTERTYPE::INDEL,"LowIndelReadPosRankSum","ReadPosRankSum","-20.0",false));
}



void
BlockerOptions::
finalize_filters() {
    if(! min_gqx.empty()) GQX_filter.reset(new FilterInfo("min-gqx",FILTERTYPE::BOTH,"LowGQX","GQX",min_gqx.c_str(),false,true,true));
}


BlockerOptions::
~BlockerOptions() {}
