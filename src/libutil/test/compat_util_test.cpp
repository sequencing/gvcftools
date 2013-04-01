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

#include "compat_util.hh"

#include <string>


BOOST_AUTO_TEST_SUITE( compat_util )

static
void
single_test_round(const double input,
                  const double expect) {

    static const double eps(0.00001);
    BOOST_CHECK_CLOSE(compat_round(input), expect, eps);
}


BOOST_AUTO_TEST_CASE( test_round ) {
    single_test_round(3.5,4.0);
    single_test_round(3.2,3.0);
    single_test_round(3.7,4.0);
    single_test_round(-1.0,-1.0);
    single_test_round(-1.2,-1.0);
    single_test_round(-1.5,-2.0);
    single_test_round(-1.7,-2.0);
}



static
void
single_test_basename(const char* input,
                     const char* expect) {

    const char* result(compat_basename(input));
    BOOST_CHECK_EQUAL(std::string(result), std::string(expect));
}


BOOST_AUTO_TEST_CASE( test_basename ) {
    single_test_basename("foo","foo");
    single_test_basename("","");

#ifdef _WIN32
    single_test_basename("\\foo","foo");
    single_test_basename("\\\\","");
    single_test_basename("\\","");
#else
    single_test_basename("/foo","foo");
    single_test_basename("//","");
    single_test_basename("/","");
#endif
}


BOOST_AUTO_TEST_SUITE_END()

