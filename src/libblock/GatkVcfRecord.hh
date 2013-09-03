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
#ifndef __GATK_VCF_RECORD_HH
#define __GATK_VCF_RECORD_HH

#include "compat_util.hh"
#include "parse_util.hh"
#include "VcfRecord.hh"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <sstream>



///
/// A value which might exist, and might be convertable to an integer
/// if it does 'convertable to an integer' in this case are ints or
/// floats which are rounded to the nearest int
///
struct MaybeInt {

    MaybeInt(const char* s)
        : IsInt((NULL != s) && ('\0' != *s) && (0 != strcmp(s,".")))
        , IntVal(0)
        , DoubleVal(0.)
    {
        if (! IsInt) return;
        if (NULL != s) StrVal = std::string(s);
        DoubleVal = parse_double(s);
        IntVal = static_cast<int>(compat_round(DoubleVal));
    }

    MaybeInt(const int i)
        : IsInt(true)
        , IntVal(i)
        , DoubleVal(i)
    {}

    bool IsNonZero() const {
        return (IsInt && (IntVal != 0));
    }

    bool IsInt;
    int IntVal;
    double DoubleVal;
    std::string StrVal;
};

std::ostream& operator<<(std::ostream& os,const MaybeInt& mi);



struct GatkVcfRecord : public VcfRecord {

    GatkVcfRecord(const istream_line_splitter& vparse)
        : VcfRecord(vparse)
        , _isgt(false)
    { }

    GatkVcfRecord(const GatkVcfRecord& orig)
        : VcfRecord(orig)
        , _isgt(orig._isgt)
        , _gt(orig._gt)
    {}

    virtual
    ~GatkVcfRecord();

    const MaybeInt& GetGQX() const {
        if (NULL == _gqx.get()) {
            MaybeInt _qual(GetQual().c_str());
            if (_qual.IsInt && GetGQ().IsInt) {
                int val = std::min(_qual.IntVal, GetGQ().IntVal);
                _gqx.reset(new MaybeInt(val));
            } else {
                _gqx.reset(new MaybeInt(""));
            }
        }
        return *_gqx;
    }

    const MaybeInt& GetGQ() const {
        if (NULL == _gq.get())
            _gq.reset(new MaybeInt(GetSampleVal("GQ")));
        return *_gq;
    }

    const MaybeInt& GetDP() const {
        if (NULL == _dp.get())
            _dp.reset(new MaybeInt(GetSampleVal("DP")));
        return *_dp;
    }

    const MaybeInt& GetMQ() const {
        if (NULL == _mq.get())
            _mq.reset(new MaybeInt(GetSampleVal("MQ")));
        return *_mq;
    }

    const std::string& GetGT() const {
        if (!_isgt) {
            GetSampleValStr("GT",_gt);
            _isgt = true;
        }
        return _gt;
    }

    bool GetIsCovered() const {
        return GetDP().IsNonZero();
    }

private:

    void
    IsSampleModified() { KillCache(); }

    void
    KillCache() {
        _gqx.reset();
        _gq.reset();
        _dp.reset();
        _mq.reset();
        _isgt=false;
    }

    mutable std::auto_ptr<MaybeInt> _gqx;
    mutable std::auto_ptr<MaybeInt> _gq;
    mutable std::auto_ptr<MaybeInt> _dp;
    mutable std::auto_ptr<MaybeInt> _mq;
    mutable bool _isgt;
    mutable std::string _gt;
};


#endif
