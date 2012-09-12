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


#include "istream_line_splitter.hh"

#include <iosfwd>


struct VcfHeaderHandler {
    VcfHeaderHandler(std::ostream& os,
                     const char* version,
                     const char* cmdline,
                     const bool is_skip_header = false)
        : _os(os)
        , _version(version)
        , _cmdline(cmdline)
        , _is_skip_header(is_skip_header)
        , _is_valid(true)
    {}

    virtual
    ~VcfHeaderHandler() {}

    bool
    process_line(const istream_line_splitter& vparse);

    void
    write_format(const char* tag,
                 const char* number,
                 const char* type,
                 const char* description) const;

protected:
    virtual
    bool
    is_skip_header_line(const istream_line_splitter& /*vparse*/) {
        return false;
    }

    virtual
    void
    process_final_header_line() {}


    std::ostream& _os;

private:
    const char* _version;
    const char* _cmdline;
    bool _is_skip_header;
    bool _is_valid;
};

#endif
