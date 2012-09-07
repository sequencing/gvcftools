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
/// an efficient (and slightly unsafe) class for basic tab-delimited files, etc...
///

/// \author Chris Saunders
///

#include "blt_exception.hh"
#include "istream_line_splitter.hh"

#include <cassert>
#include <cstring>

#include <iostream>
#include <sstream>



void
istream_line_splitter::
write_line(std::ostream& os) const {
    for(unsigned i(0);i<_n_word;++i) {
        if(i) os << _sep;
        os << word[i];
    }
    os << "\n";
}



void
istream_line_splitter::
dump(std::ostream& os) const {
    os << "\tline_no: " << _line_no << "\n";
    os << "\tline: ";
    write_line(os);
}



bool
istream_line_splitter::
parse_line() {
    _n_word=0;
    _line_no++;
    _is.getline(_buf,_buf_size);
    const unsigned len(strlen(_buf));
    if((len+1) >= _buf_size){
        std::ostringstream oss;
        oss << "ERROR: input exceeds buffer size on line_no: " << _line_no << "\n\n";
        throw blt_exception(oss.str().c_str());
    }
    
    if(! _is) {
        if(_is.eof()) { return false; } // normal eof:
        
        std::ostringstream oss;
        oss << "ERROR: Unexpected read failure in parse_line().\n";
        throw blt_exception(oss.str().c_str());
    }
    
    if(NULL == _buf) return false;
    assert(len);
    
    // do a low-level separator parse:
    {  
        char* p(_buf);
        word[0]=p;
        unsigned i(1);
        while(i<_max_word){  
            if((*p == '\n') || (*p == '\0')) break;
            if (*p == _sep) {
                *p = '\0';
                word[i++] = p+1;
            }  
            ++p;
        }
        _n_word=i;
    }
    return true;
}
