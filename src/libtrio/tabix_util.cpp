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

#include "blt_exception.hh"
#include "tabix_util.hh"

#include <sys/stat.h>

#include <iostream>
#include <sstream>


namespace {
std::ostream& log_os(std::cerr);
}



static
bool
is_prefix(const char* a,
          const char* b) {
    return (strstr(a,b)==a);
}



// this is lifted from tabix main.c:
bool
is_tabix_index(const char* f) {

    if (NULL == f) return false;

    // punt on remote case:
    if (is_prefix(f,"ftp://") || is_prefix(f,"http://")) return true;

    std::string idx(f);
    idx += ".tbi";

    struct stat stat_f,stat_idx;
    stat(f, &stat_f);
    stat(idx.c_str(), &stat_idx);
    return ( stat_f.st_mtime <= stat_idx.st_mtime );
}



void
enforce_tabix_index(const char* f) {
    if (is_tabix_index(f)) return;

    std::ostringstream oss;
    oss << "ERROR: Missing or outdated index for vcf file: " << f << "\n";
    throw blt_exception(oss.str().c_str());
}



bool
parse_tabix_region(const char* filename,
                   const char* region,
                   int& begin,
                   int& end) {

    if (NULL == region) {
        return false;
    }

    if (NULL == filename) {
        throw blt_exception("vcf filename is null ptr");
    }

    if ('\0' == *filename) {
        throw blt_exception("vcf filename is empty string");
    }

    enforce_tabix_index(filename);

    tabix_t* _tfp = ti_open(filename, 0);

    if (NULL == _tfp) {
        log_os << "ERROR: Failed to open VCF file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    // read from a specific region:
    if (ti_lazy_index_load(_tfp) < 0) {
        log_os << "ERROR: Failed to load index for vcf file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    int tid;
    const bool result(0 == ti_parse_region(_tfp->idx, region, &tid, &begin, &end));

    ti_close(_tfp);
    return result;
}



tabix_chrom_list::
tabix_chrom_list(const char* filename)
    : _nchrom(0)
    , _index(0)
    , _tptr(NULL)
    , _tfp(NULL)
{
    if (NULL == filename) {
        throw blt_exception("vcf filename is null ptr");
    }

    if ('\0' == *filename) {
        throw blt_exception("vcf filename is empty string");
    }

    tabix_t* _tfp = ti_open(filename, 0);

    if (NULL == _tfp) {
        log_os << "ERROR: Failed to open VCF file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    // read from a specific region:
    if (ti_lazy_index_load(_tfp) < 0) {
        log_os << "ERROR: Failed to load index for vcf file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    _tptr = ti_seqname(_tfp->idx,&_nchrom);
}


tabix_chrom_list::
~tabix_chrom_list() {
    if (NULL != _tptr) free(_tptr);
    if (NULL != _tfp) ti_close(_tfp);
}
