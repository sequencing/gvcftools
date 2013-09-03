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
///
/// advanced command-line option parsing for trio/twins
///
/// \author Chris Saunders
///

#include "trio_option_util.hh"

#include "boost/format.hpp"
#include "boost/program_options.hpp"


namespace po = boost::program_options;


void
validate(boost::any& v,
         const std::vector<std::string>& values,
         std::vector<info_filter>*, int) {

    if (v.empty()) {
        v = boost::any(std::vector<info_filter>());
    }
    std::vector<info_filter>* tv = boost::any_cast< std::vector<info_filter> >(&v);
    assert(NULL != tv);

    info_filter infof;

    // Extract tokens from values string vector and populate info_filter struct.
    if (values.size() != 2) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }

    infof.key = values[0];
    try {
        infof.val = boost::lexical_cast<float>(values[1]);
    } catch (const boost::bad_lexical_cast&) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }

    tv->push_back(infof);
}
