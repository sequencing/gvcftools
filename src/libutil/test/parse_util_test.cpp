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

#include "parse_util.hh"

#include <string>

BOOST_AUTO_TEST_SUITE( parse_util )


BOOST_AUTO_TEST_CASE( test_parse_int ) {
    const char* two = "2";
    const int val(parse_int(two));
    BOOST_CHECK_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str ) {
    static const char two[] = "2";
    const int val(parse_int_str(std::string(two)));
    BOOST_CHECK_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str_bad_input ) {
    static const std::string junk("ABCD");
    BOOST_CHECK_THROW(parse_int_str(junk), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()

