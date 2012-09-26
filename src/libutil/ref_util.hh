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
#ifndef __REF_UTILS_HH
#define __REF_UTILS_HH

#include "reference_contig_segment.hh"

extern "C" {
#include "faidx.h"
}

#include <string>


// specialized ref reader -- picks out many calls to individual positions:
struct samtools_char_picker {
    
    samtools_char_picker(const char* ref_file)
        : _fai(fai_load(ref_file))
    {}

    ~samtools_char_picker() { fai_destroy(_fai); }

    char
    get_char(const char* chrom,
             const int pos) const ;

private:
    faidx_t* _fai;
};


void
get_samtools_ref_seq(const char* ref_file,
                     const char* region,
                     std::string& ref_seq);

void
get_samtools_std_ref_segment(const char* ref_file,
                             const char* region,
                             reference_contig_segment& ref_seg,
                             unsigned& known_size);


struct fasta_chrom_list {

    explicit
    fasta_chrom_list(const char* filename);

    ~fasta_chrom_list();

    // returns null when names are exhuasted:
    const char* next();

private:
    int _index;
    int _nchrom;
    faidx_t* _fai;
};


#endif
