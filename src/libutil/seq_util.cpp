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
#include "seq_util.hh"

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>

namespace {
std::ostream& log_os(std::cerr);
}

void
base_error(const char* func, const char a){
    log_os << "ERROR:: Invalid base in " << func << ".\n"
           << "\t\tinvalid base (char): '" << a << "'\n"
           << "\t\tinvalid base (int): " << static_cast<int>(a) << "\n";
    exit(EXIT_FAILURE);
}



void
id_to_base_error(const unsigned i){
    log_os << "ERROR:: Invalid id in id_to_base. id: " << i << "\n";
    exit(EXIT_FAILURE);
}



bool
is_valid_seq(const char* seq){

    assert(NULL != seq);

    while(*seq !=  '\0'){
        if(! is_valid_base(*seq)) return false;
        seq++;
    }
    return true;
}



void
reverseComp(std::string& read){

    const unsigned nr(read.size());
    const unsigned nr2(nr/2);
    for(unsigned i(0);i<nr2;++i){
        const char tmp(comp_base(read[i]));
        read[i] = comp_base(read[nr-(i+1)]);
        read[nr-(i+1)] = tmp;
    }
    if(nr%2==1) read[nr2] = comp_base(read[nr2]);
}



void
get_ref_seq(const char* ref_seq_file,
            std::string& ref_seq) {

    static const unsigned buff_size(50000);
    char buff[buff_size];

    std::ifstream ref_is(ref_seq_file);

    if( ! ref_is ) {
        log_os << "ERROR:: Can't open reference sequence file: " << ref_seq_file << "\n";
        exit(EXIT_FAILURE);
    }

#ifdef SEQ_UTIL_VERBOSE
    log_os << "Reading from reference sequence file " << ref_seq_file << "\n";
#endif
    ref_is.getline(buff,buff_size);
#ifdef SEQ_UTIL_VERBOSE
    log_os << "First line: " << buff << "\n";
#endif

    ref_seq.clear();
    while(true){
        ref_is.getline(buff,buff_size);
        if(! ref_is) {
            if     (ref_is.eof()) break;
            else if(ref_is.fail()) {
                if(ref_is.bad()){
                    log_os << "ERROR:: unexpected failure while attempting to read sequence file: " << ref_seq_file << "\n";
                    exit(EXIT_FAILURE);
                }
                ref_is.clear();
            }
        }
        ref_seq += buff;
    }

#ifdef SEQ_UTIL_VERBOSE
    log_os << "Finished reading from reference sequence file " << ref_seq_file << "\n";
    log_os << "Reference sequence size: " << ref_seq.size() << "\n";
#endif
}



void
standardize_ref_seq(const char* ref_seq_file,
                    std::string& ref_seq) {

    const std::string::size_type ref_size(ref_seq.size());
    for(std::string::size_type i(0);i<ref_size;++i){
        const char old_ref(ref_seq[i]);
        char c(old_ref);
        if(islower(c)) c = toupper(c);
        if(! is_valid_base(c)){
            if(! is_iupac_base(c)) {
                log_os << "ERROR:: Unexpected character in reference sequence.\n";
                log_os << "\treference-sequence: " << ref_seq_file << "\n";
                log_os << "\tcharacter: '" << old_ref << "'\n";
                log_os << "\tcharacter-position: " << (i+1) << "\n";
                exit(EXIT_FAILURE);
            }
            c=elandize_base(c);
        }
        if(c != old_ref) ref_seq[i] = c;
    }
}
