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

#ifndef __BLOCKER_STATS_HH
#define __BLOCKER_STATS_HH


#include "stream_stat.hh"

#include <iosfwd>


struct BlockerStats {

    BlockerStats() {}

    void
    addBlock(const unsigned size,
             const double gqx_cov,
             const double dp_cov,
             const double mq_cov) {

        _block_size.add(size);

        if(size < min_cov_block()) return;
        _gqx_cov.add(gqx_cov);
        _dp_cov.add(dp_cov);
        _mq_cov.add(mq_cov);
    }

    void
    report(std::ostream& os) const;


    static
    unsigned
    min_cov_block() { return 5; }

private:
    stream_stat _block_size;

    stream_stat _gqx_cov;
    stream_stat _dp_cov;
    stream_stat _mq_cov;
};


#endif
