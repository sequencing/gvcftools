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

#ifndef __BLOCKER_VCF_HEADER_HANLDER
#define __BLOCKER_VCF_HEADER_HANLDER


#include "BlockerOptions.hh"
#include "VcfHeaderHandler.hh"

#include <string>
#include <vector>


struct BlockerVcfHeaderHandler : public VcfHeaderHandler {

    BlockerVcfHeaderHandler(const BlockerOptions& opt,
                            const char* version = NULL,
                            const char* cmdline = NULL)
        : VcfHeaderHandler(opt.outfp,version,cmdline,opt.is_skip_header)
        , _opt(opt)
    {
        static const char* rmHeaderTags[] = { "AC", "AF", "AN" };
        static const unsigned n_tags(sizeof(rmHeaderTags)/sizeof(char*));

        for(unsigned i(0);i<n_tags;++i) {
            _rmKeys.push_back(std::string("INFO=<ID=")+rmHeaderTags[i]);
        }
    }


private:
    bool
    is_skip_header_line(const istream_line_splitter& vparse);

    void
    process_final_header_line();

    void
    write_split_line(const istream_line_splitter& vparse) {
        if(_opt.is_skip_header) return;
        vparse.write_line(_os);
    }

    const BlockerOptions& _opt;
    std::vector<std::string> _rmKeys;
};

#endif
