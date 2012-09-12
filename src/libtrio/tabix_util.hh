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

#ifndef __TABIX_UTIL
#define __TABIX_UTIL

extern "C" {
#include "tabix.h"
}

// check for acceptable tabix index:
bool
is_tabix_index(const char* f);

// throw if no acceptable tabix index
void
enforce_tabix_index(const char* f);


bool
parse_tabix_region(const char* filename,
                   const char* region,
                   int& begin,
                   int& end);


struct tabix_chrom_list {

    explicit
    tabix_chrom_list(const char* filename);

    ~tabix_chrom_list();

    // returns null when names are exhuasted:
    const char* next() {
        if(_index < _nchrom) {
            return _tptr[_index++];
        } else {
            return NULL;
        }
    }

private:
    int _nchrom;
    int _index;
    const char** _tptr;
    tabix_t* _tfp;
};

#endif
