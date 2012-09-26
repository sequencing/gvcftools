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

#include "ref_util.hh"
#include "seq_util.hh"

#include <cstdlib>

#include <iostream>


namespace {
std::ostream& log_os(std::cerr);
}



char
samtools_char_picker::
get_char(const char* chrom,
         const int pos) const {

    int len; // throwaway...
    char* ref_tmp(faidx_fetch_seq(_fai,(char*)chrom,pos-1,pos-1, &len));
    if (NULL == ref_tmp) {
        log_os << "ERROR: Can't find sequence region '" << chrom << ':' << pos << '-' << pos << "' in reference file\n";
        exit(EXIT_FAILURE);
    }
    const char retval(*ref_tmp);
    free(ref_tmp);
    return retval;
}



void
get_samtools_ref_seq(const char* ref_file,
                     const char* region,
                     std::string& ref_seq) {

    faidx_t* fai(fai_load(ref_file));
    int len; // throwaway...
    char* ref_tmp(fai_fetch(fai,region, &len));
    if (NULL == ref_tmp) {
        log_os << "ERROR: Can't find sequence region '" << region << "' in reference file: '" << ref_file << "'\n";
        exit(EXIT_FAILURE);
    }
    ref_seq.assign(ref_tmp);
    free(ref_tmp);
    fai_destroy(fai);
}



void
get_samtools_std_ref_segment(const char* ref_file,
                             const char* region,
                             reference_contig_segment& ref_seg,
                             unsigned& known_size) {
    std::string& ref_seq(ref_seg.seq());
    get_samtools_ref_seq(ref_file,region,ref_seq);
    standardize_ref_seq(ref_file,ref_seq);
    known_size=0;
    const std::string::size_type ref_size(ref_seq.size());
    for(std::string::size_type i(0);i<ref_size;++i){ if(ref_seq[i]!='N') known_size++; }
}



fasta_chrom_list::
fasta_chrom_list(const char* filename)
    : _index(0)
{
    _fai=fai_load(filename);
    _nchrom=faidx_fetch_nseq(_fai);
}



fasta_chrom_list::
~fasta_chrom_list() {
    fai_destroy(_fai);
}



const char*
fasta_chrom_list::
next() {
    if(_index < _nchrom) {
        return faidx_fetch_seqname(_fai,_index++);
    } else {
        return NULL;
    }
}
