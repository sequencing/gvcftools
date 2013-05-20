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

#include "boost/test/unit_test.hpp"

#include "string_util.hh"


BOOST_AUTO_TEST_SUITE( string_util )


BOOST_AUTO_TEST_CASE( test_split_string ) {
    static const std::string test("A:B:C");
    static const char delim(':');

    std::vector<std::string> words;
    words.push_back("DECOY");
    split_string(test,delim,words);
    BOOST_CHECK_EQUAL(words.size(), 3);
    BOOST_CHECK_EQUAL(words[0], "A");
}

BOOST_AUTO_TEST_CASE( test_split_string_cstr ) {
    static const char* test("A:B:C");
    static const char delim(':');

    std::vector<std::string> words;
    words.push_back("DECOY");
    split_string(test,delim,words);
    BOOST_CHECK_EQUAL(words.size(), 3);
    BOOST_CHECK_EQUAL(words[0], "A");
}

BOOST_AUTO_TEST_SUITE_END()

