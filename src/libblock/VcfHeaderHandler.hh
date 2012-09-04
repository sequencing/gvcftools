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

#ifndef __VCF_HEADER_HANLDER
#define __VCF_HEADER_HANLDER


#include "BlockerOptions.hh"
#include "istream_line_splitter.hh"

#include <iosfwd>
#include <string>


struct VcfHeaderHandler {

    VcfHeaderHandler(const BlockerOptions& opt,
                     const char* version = NULL,
                     const char* cmdline = NULL)
        : _opt(opt)
        , _version(version)
        , _cmdline(cmdline)
        , _outfp(opt.outfp)
        , _is_valid(true)
    {
        static const char* rmHeaderTags[] = { "AC", "AF", "AN" };
        static const unsigned n_tags(sizeof(rmHeaderTags)/sizeof(char*));

        for(unsigned i(0);i<n_tags;++i) {
            _rmKeys.push_back(std::string("INFO=<ID=")+rmHeaderTags[i]);
        }
    }

    bool
    process_line(const istream_line_splitter& vparse);

private:

    void
    write_split_line(const istream_line_splitter& vparse);


    const BlockerOptions& _opt;
    const char* _version;
    const char* _cmdline;
    std::ostream& _outfp;
    bool _is_valid;

    std::vector<std::string> _rmKeys;
};

#endif
