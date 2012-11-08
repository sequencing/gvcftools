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
#ifndef __VCF_RECORD_BLOCKER_HH
#define __VCF_RECORD_BLOCKER_HH

#include "BlockerOptions.hh"
#include "BlockVcfRecord.hh"

#include <string>



/// Accumulate vcf_record objects and determine when these will be combined into block records
///
struct VcfRecordBlocker {

    VcfRecordBlocker(const BlockerOptions& opt)
        : _opt(opt)
        , _blockCvcfr(opt,_stats)
        , _is_highDepth(false)
        , _bufferStartPos(0)
        , _bufferEndPos(0)
        , _lastNonindelPos(0)
    {}

    /// Process and print any remaining blocks
    ~VcfRecordBlocker();

    /// Submit next vcf record for printing or blocking
    ///
    void Append(GatkVcfRecord& record)
    {
        if(IsSkipRecord(record)) return;

        GroomInputRecord(record);
        AccumulateRecords(record);
    }

private:

    void WriteBlockCvcfr() {
        _blockCvcfr.Write(_opt.outfp);
        _blockCvcfr.Reset();
    }

    // write a single non-blockable record:
    void WriteThisCvcfr(GatkVcfRecord& record) {
        if (record.GetGQX().IsInt) {
            record.SetSampleVal("GQX",_intstr.get32(record.GetGQX().IntVal));
        }

        record.WriteUnaltered(_opt.outfp);
    }

    bool IsRecordInCurrentBlock(GatkVcfRecord& record) {
        return _blockCvcfr.Test(record);
    }

    void JoinRecordToBlock(GatkVcfRecord& record){
        _blockCvcfr.Add(record);
    }

    bool
    IsSkipRecord(GatkVcfRecord& record) {
        // check to make sure this isn't a reference block record from GATK -- I believe these are created when
        // an indel is evaluated but the reference allele is chosen, it would be nice to incorporate these
        // into the final gVCF but until there's a policy it's cleaner to filter such cases:
        if(record.IsNonvariantBlock()) return true;

        // besides the non-variant blocks, GATK will infrequently output the same reference site twice, so filter
        // thest cases out if found:
        if(! record.IsIndel()) {
            const unsigned pos(record.GetPos());
            if(pos <= _lastNonindelPos) return true;
            _lastNonindelPos = pos;
        }


        return false;
    }

    void GroomInputRecord(GatkVcfRecord& record);

    // accumulate all contiguous regions where sites or indels overlap with other indels:
    //
    void AccumulateRecords(GatkVcfRecord& record) {

        const bool is_indel(record.IsIndel());
        const int pos(record.GetPos());

        if(is_indel) {
            // update start and end pos:
            const int startPos(pos+1);
            const int endPos(pos+static_cast<int>(record.GetRef().size())-1);

            // We use an expanded definition of indel overlap between
            // indels here such that an overlapping insertion,
            // insertion, deletion (in that order) would all fall into
            // the same buffer group.
            
            const bool is_in_indel(! (startPos+1>_bufferEndPos || endPos+1<_bufferStartPos));
            //            const bool is_in_indel((pos>=_bufferStartPos) && (pos<=_bufferEndPos));
            if(is_in_indel) {
                _bufferStartPos=std::min(_bufferStartPos,startPos);
            } else {
                _bufferStartPos=startPos;
            }

            _bufferEndPos=std::max(_bufferEndPos,endPos);
            
            if(!(_recordBuffer.empty() || is_in_indel)){
                ProcessRecordBuffer();
            }
            _indelIndex.push_back(_recordBuffer.size());
            _recordBuffer.push_back(record);
        } else {
            const bool is_in_indel((pos>=_bufferStartPos) && (pos<=_bufferEndPos));
            if(is_in_indel) {
                _recordBuffer.push_back(record);
            } else {
                if(!_recordBuffer.empty()) {
                    ProcessRecordBuffer();
                }
                ProcessRecord(record);
            }
        }
    }

    void GroomRecordBuffer();

    // make changes to records in the record buffer according to
    // overlapping indel/site handling policies (via
    // GroomRecordBuffer), then submit all buffered records for block
    // compression/write:
    void ProcessRecordBuffer();


    // After all input modifications are finished, send record on to
    // be joined to a block or printed:
    void ProcessRecord(GatkVcfRecord& record) {

        if (!IsVcfRecordBlockable(record)) {
            WriteBlockCvcfr();
            WriteThisCvcfr(record);
            return;
        }

        if (!IsRecordInCurrentBlock(record)) {
            WriteBlockCvcfr();
        }
        JoinRecordToBlock(record);
    }

    static
    void
    AddFilterSet(GatkVcfRecord& record,
                 const std::vector<FilterInfo>& filters) {

        const bool is_indel(record.IsIndel());
        const unsigned fs(filters.size());
        for(unsigned i(0);i<fs;++i) {
            const FILTERTYPE::index_t ft(filters[i].filter_type);
            if       (ft == FILTERTYPE::SITE) {
                if(is_indel) continue;
            } else if(ft == FILTERTYPE::INDEL) {
                if(! is_indel) continue;
            }
            AddFilter(record,filters[i]);
        }
    }

    static
    void
    AddFilter(GatkVcfRecord& record,
              const FilterInfo& filter);

    /// Certain vcf records can *never* be compressed -- such as variants
    /// and annotated sites
    bool
    IsVcfRecordBlockable(GatkVcfRecord& record) {
        if (record.GetId() != ".") return false;
        if (record.IsVariant()) return false;
        if (record.GetRef().size() != 1) return false;
        //if (record.GetInfo().Count != 0)
        //    return false; // might have to take this one out eventually

        const std::string& gt = record.GetGT();
        if ((!gt.empty()) && (gt != "./.") && (gt != ".") && (gt != "0/0") && (gt != "0")) return false;

        // AD from GATK uses unfiltered counts, for this reason we use
        // info DP (unfiltered) instead of sample DP (filtered)
        const MaybeInt info_dp(record.GetInfoVal("DP"));
        const MaybeInt ad(record.GetSampleVal("AD"));
        if(ad.IsInt && info_dp.IsInt) {
            const double reffrac(static_cast<double>(ad.IntVal)/static_cast<double>(info_dp.IntVal));
            if((reffrac+_opt.min_nonref_blockable.numval()) <= 1.0) return false;
        }

        return true;
    }

private:
    const BlockerOptions& _opt;
    BlockVcfRecord _blockCvcfr;

    std::string _lastChrom;
    bool _is_highDepth;
    double _highDepth;

    int _bufferStartPos,_bufferEndPos; // buffer all records on [Start,End]
    std::vector<GatkVcfRecord> _recordBuffer; // buffer positions crossed by deletions or other indel events
    std::vector<unsigned> _indelIndex; // record index of records in buffer which are indels

    unsigned _lastNonindelPos;

    //tmp catch for gt parsing:
    std::vector<int> _gti;
    //obj for fast int->str
    stringer<int> _intstr;

    BlockerStats _stats;
};


#endif
