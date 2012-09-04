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
#ifndef __REFERENCE_CONTIG_SEGMENT_HH
#define __REFERENCE_CONTIG_SEGMENT_HH

#include "pos_type.hh"

#include <string>




/// Manages a partial reference sequence segment
///
/// This object holds the reference sequence specified by the current
/// runs begin and end range, plus some padding on each side. To get 
/// this integrated into the current code as quickly as possible it
/// currently exposes the internal string object holding the sequence
/// data. When time allows this will be restricted so that a compressed
/// internal object can be used.
///
struct reference_contig_segment {

    reference_contig_segment()
        : _offset(0)
    {}

    char
    get_base(const pos_t pos) const {
        if(pos<_offset || pos>=end()) return 'N';
        return _seq[pos-_offset];
    }

    std::string& seq() { return _seq; }
    const std::string& seq() const { return _seq; }

    void
    set_offset(const pos_t offset) {
        _offset=offset;
    }

    pos_t
    end() const { return _offset+_seq.size(); }

private:

    pos_t _offset;
    std::string _seq;
};


#endif
