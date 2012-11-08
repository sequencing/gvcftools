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
#ifndef __BLOCK_VCF_RECORD_HH
#define __BLOCK_VCF_RECORD_HH


#include "BlockerOptions.hh"
#include "BlockerStats.hh"
#include "GatkVcfRecord.hh"
#include "stream_stat.hh"
#include "stringer.hh"

#include <cassert>
#include <cstdio>

#include <memory>


/// stores vcf records representing contiguous blocks of non-variant sites
///
struct BlockVcfRecord {

    BlockVcfRecord(const BlockerOptions& opt,
                   BlockerStats& stats)
        : _opt(opt)
        , _fracTol(opt.nvopt.BlockFracTol.numval())
        , _absTol(opt.nvopt.BlockAbsTol)
        , _count(0)
        , _stats(stats)
    {}

    ~BlockVcfRecord();

    void Reset() {
        _baseCvcfr.reset();
        _count=0;
        _blockGQX.reset();
        _blockDP.reset();
        _blockMQ.reset();
    }

    /// determine if new record can be incorporated into the current block
    bool Test(GatkVcfRecord& cvcfr) const {

        if (_count == 0) return true;

        // check if chrom matches and pos is +1 from end record:
        if (cvcfr.GetChrom() != _baseCvcfr->GetChrom())
            return false;
        if (cvcfr.GetPos() != (_baseCvcfr->GetPos() + _count))
            return false;

        // does the filter field match?
        if (cvcfr.GetFilter() != _baseCvcfr->GetFilter())
            return false;

        // does the gt field match?
        if (cvcfr.GetGT() != _baseCvcfr->GetGT())
            return false;

        // special check for no-coverage regions:
        if (_baseCvcfr->GetIsCovered() != cvcfr.GetIsCovered())
            return false;

        // none of the checks below apply to no-coverage regions
        if (!_baseCvcfr->GetIsCovered())
            return true;

        // test gq
        if (!IsNewValueBlockable(cvcfr.GetGQX(),_baseCvcfr->GetGQX(),
                                 _blockGQX,_fracTol,_absTol))
                return false;

        if (!IsNewValueBlockable(cvcfr.GetDP(),_baseCvcfr->GetDP(),
                                 _blockDP,_fracTol,_absTol))
                return false;

        if (!IsNewValueBlockable(cvcfr.GetMQ(),_baseCvcfr->GetMQ(),
                                 _blockMQ,_fracTol,_absTol))
                return false;

        return true;
    }

    void
    Add(GatkVcfRecord& cvcfr) {
        if (_count == 0)
            _baseCvcfr.reset(new GatkVcfRecord(cvcfr));

        if (cvcfr.GetGQX().IsInt)
            _blockGQX.add(cvcfr.GetGQX().IntVal);
        if (cvcfr.GetDP().IsInt)
            _blockDP.add(cvcfr.GetDP().IntVal);
        if (cvcfr.GetMQ().IsInt)
            _blockMQ.add(cvcfr.GetMQ().IntVal);
        
        _count += 1;
    }

    void
    Write(std::ostream& os){

        if (_count == 0) return;

        // store preserved sample information:
        const std::string& gt = _baseCvcfr->GetGT();

        // get this value before clearing record:
        const bool is_covered(_baseCvcfr->GetIsCovered());
        // clear most of record:
        _baseCvcfr->ClearInfo();
        _baseCvcfr->ClearSample();
        _baseCvcfr->SetQual(NULL);

        // add some fields back:
        _baseCvcfr->SetSampleVal("GT", gt.c_str());

        if (_count > 1) {
            static const unsigned buff_size(32);
            char buff[buff_size];

            const int end(_baseCvcfr->GetPos() + _count - 1);
            const int write_size(snprintf(buff,buff_size,"END=%i",end));
            assert((write_size>=0) && (write_size < static_cast<int>(buff_size)));
            _baseCvcfr->AppendInfo(buff);
        }

        // insert averaged values for covered blocks:
        if (is_covered) {
            bool isAvg(false);

            UpdateBlock("DP",_blockDP,isAvg);
            UpdateBlock("GQX",_blockGQX,isAvg);
            UpdateBlock("MQ",_blockMQ,isAvg);

            if (isAvg) {
                const std::string& label(_opt.nvopt.BlockavgLabel);
                _baseCvcfr->AppendInfo(label.c_str());
            }
        }

        _stats.addBlock(_count);

        _baseCvcfr->WriteUnaltered(os);
    }


private:

    void
    UpdateBlock(const char* label,
                stream_stat& block,
                bool& isAvg) {

        static const char* unknown = ".";
        const char* printptr(unknown);

        if (! block.empty()) {
            if (block.size() > 1) isAvg = true;
            const int min(static_cast<int>(compat_round(block.min())));
            printptr=_intstr.get32(min);
        }
        _baseCvcfr->SetSampleVal(label,printptr);
    }
    
    static
    bool 
    IsNewValueBlockable(const MaybeInt& newval,
                        const MaybeInt& oldval,
                        const stream_stat& ss,
                        const double fracTol,
                        const int absTol) {
        if (!(newval.IsInt && oldval.IsInt)) {
                return (newval.StrVal == oldval.StrVal);
        }
        return IsNewValueBlockableInternal(newval.IntVal, ss, fracTol, absTol);
    }

    // should be called only after new/old null state has been queried
    static
    bool
    IsNewValueBlockableInternal(const int newval,
                                const stream_stat& ss,
                                const double fracTol,
                                const int absTol) {

        stream_stat ss2(ss);
        ss2.add(newval);
        return CheckBlockTolerance(ss2, fracTol, absTol);
    }

    static
    bool
    CheckBlockSingleTolerance(const stream_stat& ss,
                              const int min,
                              const int tol) {
        return ((min + tol) >= ss.max());
    }

    /// <summary>
    /// check that int(block_mean) is with either abs_tol of min/max or mean*frac_tol of min/max
    /// </summary>
    static 
    bool
    CheckBlockTolerance(const stream_stat& ss,
                        const double fracTol,
                        const int absTol) {

        const int min(static_cast<int>(compat_round(ss.min())));
        if (CheckBlockSingleTolerance(ss, min, absTol)) return true;
        const int ftol(static_cast<int>(std::floor(min * fracTol)));
        if (ftol <= absTol) return false;
        return CheckBlockSingleTolerance(ss, min, ftol);
    }

    const BlockerOptions& _opt;
    const double _fracTol;
    const int _absTol;
    std::auto_ptr<GatkVcfRecord> _baseCvcfr;
    int _count;
    BlockerStats _stats;

    stream_stat _blockGQX;
    stream_stat _blockDP;
    stream_stat _blockMQ;

    // fast int->str util:
    stringer<int> _intstr;
};





#endif
