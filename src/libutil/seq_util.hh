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
#ifndef __SEQ_UTIL_HH
#define __SEQ_UTIL_HH

#include <string>


enum { N_BASE=4 };

void
base_error(const char* func,
           const char a);

inline
unsigned
base_to_id(const char a) {
    switch (a) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:
        base_error("base_to_id",a);
        return 4;
    }
}

void
id_to_base_error(const unsigned i);

inline
char
id_to_base(const unsigned i) {
    static const char base[] = "ACGT";

    if (i>=N_BASE) id_to_base_error(i);
    return base[i];
}



/// valid in the ELAND sense [ACGTN]
inline
bool
is_valid_base(char a) {
    switch (a) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
    case 'N': return true;
    default : return false;
    }
}

inline
bool
is_iupac_base(char a) {
    switch (a) {
    case 'A':
    case 'C':
    case 'G':
    case 'U':
    case 'T':
    case 'R':
    case 'Y':
    case 'S':
    case 'W':
    case 'K':
    case 'M':
    case 'B':
    case 'D':
    case 'H':
    case 'V':
    case '.':
    case '-':
    case 'N': return true;
    default : return false;
    }
}

/// valid in the ELAND sense [ACGTN]
bool
is_valid_seq(const char* seq);

inline
char
elandize_base(char a) {
    switch (a) {
    case 'A': return 'A';
    case 'C': return 'C';
    case 'G': return 'G';
    case 'U':
    case 'T': return 'T';
    case 'R':
    case 'Y':
    case 'S':
    case 'W':
    case 'K':
    case 'M':
    case 'B':
    case 'D':
    case 'H':
    case 'V':
    case '.':
    case '-':
    case 'N': return 'N';
    default:
        base_error("elandize_base",a);
        return 'N';
    }
}

inline
char
comp_base(char a) {
    switch (a) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    default:
        base_error("comp_base",a);
        return 'N';
    }
}

void
reverseComp(std::string& read);

void
get_ref_seq(const char* ref_seq_file,
            std::string& ref_seq);

/// Standardize reference sequence to [ACGTN]. Fail when non-IUPAC
/// character is found.
void
standardize_ref_seq(const char* ref_seq_file,
                    std::string& ref_seq);

#endif
