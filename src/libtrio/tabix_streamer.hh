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
/// No frills object wrapper for tabix region or whole file line iteration
///

/// \author Chris Saunders
///
#ifndef __TABIX_STREAMER_HH
#define __TABIX_STREAMER_HH

#include "blt_exception.hh"

extern "C" {
#include "tabix.h"
}

#include <iosfwd>
#include <string>


struct tabix_streamer {

    // optionally provide a BAM header to validatr vcf chromosome names against
    //
    explicit
    tabix_streamer(const char* filename,
                   const char* region = NULL);
    
    ~tabix_streamer();

    bool next();

    // return NULL on finish
    char* getline() const {
        if(_is_record_set) return _linebuf;
        else               return NULL;
    }

    const char* name() const { return _stream_name.c_str(); }

    unsigned record_no() const { return _record_no; }

    void report_state(std::ostream& os) const;

private:
    bool _is_record_set;
    bool _is_stream_end;
    unsigned _record_no;
    std::string _stream_name;

    tabix_t* _tfp;
    ti_iter_t _titer;

    char* _linebuf;
};


#endif
