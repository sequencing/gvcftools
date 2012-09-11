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

#include "VcfHeaderHandler.hh"

#include <cstring>

#include <iostream>


bool
VcfHeaderHandler::
process_line(const istream_line_splitter& vparse){
    if(! _is_valid) return false;
    
    const unsigned nw(vparse.n_word());
    if ((nw == 0) || (vparse.word[0][0] != '#')) {
        // shouldn't get here unless the header is missing, etc...
        _is_valid = false;
        return _is_valid;
    }
    
    const bool is_last(0 == strcmp(vparse.word[0],"#CHROM"));
    if(is_last) _is_valid = false;
    
    if(_is_skip_header) return true;
    
    if(! is_last) {
        if(is_skip_header_line(vparse)) return _is_valid;
    } else {
        if(NULL != _version) {
            _os << "##gvcftools_version=\"" << _version << "\"\n";
        }
        if(NULL != _cmdline) {
            _os << "##gvcftools_cmdline=\"" << _cmdline << "\"\n";
        }
        
        process_final_header_line();
    }

    vparse.write_line(_os);
    return true;
}



void
VcfHeaderHandler::
write_format(const char* tag,
             const char* number,
             const char* type,
             const char* description) const {

    _os << "##FORMAT=<ID=" << tag
        << ",Number=" << number
        << ",Type=" << type
        << ",Description=\"" << description << "\">\n";
}
