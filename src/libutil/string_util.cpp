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


#include "string_util.hh"

#include <cstring>

#include <iostream>



void
split_string(const char* str,
             const char delimiter,
             std::vector<std::string>& v) {

    while(true) {
        const char* next(strchr(str,delimiter));
        if((NULL == next) || (delimiter == '\0')) {
            v.push_back(std::string(str));
            return;
        }
        v.push_back(std::string(str,next-str));
        str = next+1;
    }
}



void
split_string(const std::string& str,
             const char delimiter,
             std::vector<std::string>& v) {

    size_t start(0);
    while(true) {
        size_t next(str.find(delimiter,start));
        v.push_back(std::string(str.substr(start,next-start)));
        if(next == std::string::npos) return;
        start = next+1;
    }
}



bool
split_match(const std::string& str,
            const char delimiter,
            const char* needle) {

    size_t start(0);
    while(true) {
        size_t next(str.find(delimiter,start));
        if(0 == str.compare(start,next-start,needle)) return true;
        if(next == std::string::npos) break;
        start = next+1;
    }
    return false;
}



static
void
dump_result(const char* str,
            const std::vector<std::string>& v) {
    std::ostream& os(std::cerr);

    os << "input test string '" << str << "'\n";
    for(unsigned i(0);i<v.size();++i) {
        os << "word: '" << v[i] << "'\n";
    }
    os << "\n";
}


static
void
test_split(const char* test) {

    std::vector<std::string> v1;
    split_string(test,',',v1);
    dump_result(test,v1);

    std::vector<std::string> v2;
    std::string s(test);    
    split_string(s,',',v2);
    dump_result(test,v2);
}



static
void
test_match(const char* str,
           const char delimiter,
           const char* needle) {
    std::ostream& os(std::cerr);

    std::string s(str);
    const bool res = split_match(s,delimiter,needle);
    os << "input test string '" << str << "'\n";
    os << "find: '" << needle << "'\n";
    os << "result: " << res << "\n";
    os << "\n";
}



void
run_string_util_tests() {
    test_split(",1,3,5");
    test_split("");
    test_split(",");
    test_split("1234");

    test_match("foo,boo,",',',"foo");
    test_match("foo,boo,",',',"xoo");
}


//int main() { run_string_util_tests(); }

