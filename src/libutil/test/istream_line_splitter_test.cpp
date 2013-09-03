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

#include "istream_line_splitter.hh"

#include <sstream>
#include <string>


BOOST_AUTO_TEST_SUITE( istream_line_splitter_test )


BOOST_AUTO_TEST_CASE( itest_istream_line_splitter )
{

    std::string test_input("1\t2\t3\t4\n11\t22\t33\t44\n");
    std::istringstream iss(test_input);

    istream_line_splitter dparse(iss);

    int line_no(0);
    while (dparse.parse_line()) {
        line_no++;
        static const unsigned expected_col_count(4);
        BOOST_CHECK_EQUAL(dparse.n_word(),expected_col_count);
        if       (1==line_no) {
            BOOST_CHECK_EQUAL(std::string(dparse.word[1]),std::string("2"));
        } else if (2==line_no) {
            BOOST_CHECK_EQUAL(std::string(dparse.word[1]),std::string("22"));
        }
    }
}


static
void
check_long_line(const int init_buffer_size) {

    std::string test_input("1ABCDEFGHIJKLMNOPQRSTUVWXYZ\t2\t3\t4ABCDEFG\n11\t22\t33\t44XYZ\n");
    std::istringstream iss(test_input);

    istream_line_splitter dparse(iss,init_buffer_size);

    int line_no(0);
    while (dparse.parse_line()) {
        line_no++;
        static const unsigned expected_col_count(4);
        BOOST_CHECK_EQUAL(dparse.n_word(),expected_col_count);
        if       (1==line_no) {
            BOOST_CHECK_EQUAL(std::string(dparse.word[0]),std::string("1ABCDEFGHIJKLMNOPQRSTUVWXYZ"));
            BOOST_CHECK_EQUAL(std::string(dparse.word[3]),std::string("4ABCDEFG"));
        } else if (2==line_no) {
            BOOST_CHECK_EQUAL(std::string(dparse.word[2]),std::string("33"));
        }
    }
}


BOOST_AUTO_TEST_CASE( itest_istream_line_splitter_long_line )
{
    check_long_line(2);
    check_long_line(41);
}

BOOST_AUTO_TEST_SUITE_END()

