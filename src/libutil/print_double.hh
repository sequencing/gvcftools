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
#ifndef __PRINT_DOUBLE_HH
#define __PRINT_DOUBLE_HH

#include "parse_util.hh"

#include "boost/program_options.hpp"

#include <iosfwd>
#include <string>
#include <vector>


/// holds floating-point values with exact control of pretty print:
///
struct print_double {
    print_double(const char* init = NULL)
        : _numval(0)
    {
        if(NULL == init) {
            _strval="0";
            return;
        }
        update(init);
    }

    void
    update(const char* str) {
        _strval=str;
        _numval=parse_double_str(_strval);
    }

    const std::string&
    strval() const { return _strval; }

    double
    numval() const { return _numval; }

private:
    std::string _strval;
    double _numval;
};


std::ostream& operator<<(std::ostream& os,const print_double& pd);


// this allows print_double to work with boost program_options
void validate(boost::any& v,
              const std::vector<std::string>& values,
              print_double*, int);

#endif
