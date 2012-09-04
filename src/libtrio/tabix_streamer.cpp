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

#include "tabix_streamer.hh"

#include <cassert>
#include <cstdlib>
//#include <sys/stat.h>

#include <iostream>
#include <set>
#include <string>

namespace {
std::ostream& log_os(std::cerr);
}


tabix_streamer::
tabix_streamer(const char* filename,
               const char* region)
    : _is_record_set(false)
    , _is_stream_end(false)
    , _record_no(0)
    , _stream_name(filename)
    , _tfp(NULL)
    , _titer(NULL)
    , _linebuf(NULL)
{

    if(NULL == filename){
        throw blt_exception("vcf filename is null ptr");
    }

    if('\0' == *filename){
        throw blt_exception("vcf filename is empty string");
    }

    _tfp = ti_open(filename, 0);

    if(NULL == _tfp) {
        log_os << "ERROR: Failed to open VCF file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    // read from a specific region:
    if (ti_lazy_index_load(_tfp) < 0) {
        log_os << "ERROR: Failed to load index for vcf file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    if(NULL == region) {
        // read the whole VCF file:
        _titer = ti_query(_tfp, 0, 0, 0);
        return;
    }

    int tid,beg,end;
    if (ti_parse_region(_tfp->idx, region, &tid, &beg, &end) == 0) {
        _titer = ti_queryi(_tfp, tid, beg, end);
    } else {
        _is_stream_end=true;
    }
}



tabix_streamer::
~tabix_streamer() {
    if(NULL != _titer) ti_iter_destroy(_titer);
    if(NULL != _tfp) ti_close(_tfp);
}



bool
tabix_streamer::
next() {
    if(_is_stream_end || (NULL==_tfp) || (NULL==_titer)) return false;

    int len;
    _linebuf = (char*) ti_read(_tfp, _titer, &len);
    
    _is_stream_end=(NULL == _linebuf);
    _is_record_set=(! _is_stream_end);
    if(_is_record_set) _record_no++;
    
    return _is_record_set;
}



void
tabix_streamer::
report_state(std::ostream& os) const {

    const char* line(getline());
    os << "\tvcf_stream_label: " << name() << "\n";
    if(NULL != line){
        os << "\tvcf_stream_record_no: " << record_no() << "\n"
           << "\tvcf_record: '" << line << "'\n";
    } else {
       os << "\tno vcf record currently set\n";
    }
}
