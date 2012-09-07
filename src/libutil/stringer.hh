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

#ifndef __STRINGER_HH
#define __STRINGER_HH

#include "compat_util.hh"

#include <cstdio>


/// String conversion utility which is harder-to-use but faster than stringstream/lexical_cast
///
/// Safety notes:
/// 1) client must create one object for each thread
/// 2) The string pointer returned will be invalid at the next conversion call to stringer
///
struct stringer {

    const char*
    itos_32(const int i) {
        static const unsigned buff_size(32);
        const int write_size(snprintf(_buff32,buff_size,"%i",i));
        if((write_size<0) || (write_size >= static_cast<int>(buff_size))) {
            itos_32_error(i,write_size);
        }
        return _buff32;
    }

private:
    void itos_32_error(const int i, const int write_size);

    char _buff32[32];
};


#endif
