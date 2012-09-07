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


#include "parse_util.hh"
#include "VcfRecordBlocker.hh"

#include "boost/lexical_cast.hpp"

//#define VDEBUG

#ifdef VDEBUG 
#include <iostream>
#endif



static
bool
checked_double_parse(const char* s,
                     double& val) {
    if((NULL != s) && ('\0' != *s) && (0 != strcmp(s,"."))) return false;
    val=parse_double(s);
    return true;
}



void
VcfRecordBlocker::
GroomInputRecord(GatkVcfRecord& record) {
    // handle special filters:

    // GQX needs to be handled separately because it is derived,
    // rather than actually in, the input vcf record:
    if(NULL != _opt.GQX_filter.get()) {
        const MaybeInt& gqx(record.GetGQX());
        if ((!gqx.IsInt) || (gqx.DoubleVal < _opt.GQX_filter->thresh.numval)) {
            record.AppendFilter(_opt.GQX_filter->label.c_str());
        }
    }
    
    // high depth filter:
    if(_opt.is_chrom_depth()) {
        // filter for high depth:
        const std::string& thisChrom(record.GetChrom());
        if ((_lastChrom.empty()) || (_lastChrom != thisChrom)) {
            
            _is_highDepth=(0 != _opt.ChromDepth.count(thisChrom));
            if(_is_highDepth) {
                _highDepth = _opt.ChromDepth.find(thisChrom)->second * _opt.max_chrom_depth_filter_factor.numval;
            }
            _lastChrom = thisChrom;
        }
        
        if (_is_highDepth) {
            const char* dp(record.GetSampleVal("DP"));
            if ((NULL != dp) && (parse_double(dp) > _highDepth)){
                record.AppendFilter(_opt.max_chrom_depth_filter_tag.c_str());
            }
        }
    }

    // handle all other filters:
    AddFilterSet(record,_opt.filters);

    // handle newer GATK-input case where "." is used for filter field instead of "PASS"
    if (record.GetFilter().empty()) { record.PassFilter(); }

    // remove standard pop-gen info tags from variant and non-variant records:
    record.DeleteInfoKeyVal("AC");
    record.DeleteInfoKeyVal("AF");
    record.DeleteInfoKeyVal("AN");

    // transfer MQ over to a sample value for block averaging
    MaybeInt mqVal(record.GetInfoVal("MQ"));
    if (! mqVal.StrVal.empty()) {
        if (mqVal.IsInt) {
            static const unsigned buff_size(32);
            char buff[buff_size];
            const int write_size(snprintf(buff,buff_size,"%i",mqVal.IntVal));
            assert((write_size>=0) && (write_size < static_cast<int>(buff_size)));

            record.SetSampleVal("MQ", buff);
        } else {
            record.SetSampleVal("MQ", mqVal.StrVal.c_str());
        }
    }
}



void
VcfRecordBlocker::
AddFilter(GatkVcfRecord& record,
          const FilterInfo& filter) {

    bool is_filter(false);
    const char* tagval(NULL);
    if(filter.is_sample_value) {
        tagval = record.GetSampleVal(filter.tag.c_str());
    } else {
        tagval = record.GetInfoVal(filter.tag.c_str());
    }
    const MaybeInt val(tagval);
    if (!val.IsInt) {
        if(filter.is_filter_if_missing) {
            is_filter=true;
        }
    } else {
        if ((filter.is_max_thresh && (val.DoubleVal > filter.thresh.numval)) ||
            ((!filter.is_max_thresh) && (val.DoubleVal < filter.thresh.numval)))
            is_filter=true;
    }
    if(is_filter) record.AppendFilter(filter.label.c_str());
}



struct double_info {

    double_info()
        : is_valid(false)
        , val(0)
        , str(NULL)
    {}

    bool is_valid;
    double val;
    const char* str;
};



// info recorded for each indel spanning region:
struct region_info {

    region_info()
        : copyn(0)
    {}

    // record information on our confidence in the indel call. The
    // indel confidence represents an upper-bound on our confidence of
    // site calls within the indel.
    std::vector<std::string> filters;
    double_info qual;
    double_info gq;

    // For the whole region, record whether we expect this to be
    // hemizgous (a heterozygous indel) or not-covered (a homozygous
    // indel). For now, overlapping indels are treated the same as
    // not-covered, because we assume that alignment is not specified
    // in the input VCF.
    //
    // TODO: If overlapping indels are supplied with CIGAR alignments
    // for each alternate allele, then it's possible to make a site
    // map of hemizygous sites for these cases.
    //const unsigned region_size(_bufferEndPos+1-static_cast<int>(_bufferStartPos));
    unsigned copyn;
};



// function used for sites inside of homozygous deletions and 
// unresolvable sites inside of heterozygous deletions:
static
void
set_record_to_unknown_gt(GatkVcfRecord& record) {
    record.SetQual(".");
    record.DeleteSampleKeyVal("PL");
    record.DeleteSampleKeyVal("GQ");
    record.DeleteSampleKeyVal("GQX");
    record.SetSampleVal("GT",".");
}



typedef std::pair<unsigned,char> refedit;


static
void
adjust_overlap_record(const BlockerOptions& opt,
                      const region_info& rinfo,
                      const unsigned,
                      GatkVcfRecord& record,
                      bool&,
                      std::vector<refedit>&) {

    // apply filters:
    const unsigned n_filt(rinfo.filters.size());
    for(unsigned filt_index(0);filt_index<n_filt;++filt_index) {
        record.AppendFilter(rinfo.filters[filt_index].c_str());
    }
    
    //apply quality minimums:
    if(rinfo.qual.is_valid) {
        double record_qual(0.);
        const bool is_valid(checked_double_parse(record.GetQual().c_str(),record_qual));
        if(is_valid && (rinfo.qual.val<record_qual)) {
            record.SetQual(rinfo.qual.str);
        }
    }
    
    if(rinfo.gq.is_valid) {
        double record_gq(0.);
        const bool is_valid(checked_double_parse(record.GetSampleVal("GQ"),record_gq));
        if(is_valid && (rinfo.gq.val<record_gq)) {
            record.SetSampleVal("GQ",rinfo.gq.str);
        }
    }
    
    // change gt conflict status based on region_copyn
    assert(rinfo.copyn<2);

    if(rinfo.copyn==1) {
        std::vector<int> gti;
        if(! record.GetGT().empty()) {
            parse_gt(record.GetGT().c_str(),gti);
        }

        if(gti.size() == 2) {
            if(gti[0]==gti[1]) {
                if       (gti[0]>=0) {
                    std::ostringstream oss;
                    oss << gti[0];
                    record.SetSampleVal("GT",oss.str().c_str());
                    record.DeleteSampleKeyVal("PL");
                } else {
                    set_record_to_unknown_gt(record);
                }
            } else {
                set_record_to_unknown_gt(record);
                record.AppendFilter(opt.site_conflict_label.c_str());
            }
        } else if(gti.size() != 1) {            
            set_record_to_unknown_gt(record);
        }
    } else {
        set_record_to_unknown_gt(record);
    }
}



// modify overlapping site and indel records to be self-consistent:
void
VcfRecordBlocker::
GroomRecordBuffer() {

    const unsigned n_records(_recordBuffer.size());

#ifdef VDEBUG 
    if(true) {
        std::cerr << "VDEBUG input: indel count: " << _indelIndex.size() << "\n";
        for(unsigned i(0);i<n_records;++i) {
            _recordBuffer[i].WriteUnaltered(std::cerr);
        }
    }
#endif

    // create a map of 'covered' ploidy through the indel region based
    // on the first indel, any additional inside of the first must be
    // conflict:
    region_info rinfo;

    if(_indelIndex.size() > 1) {
        rinfo.filters.push_back(_opt.indel_conflict_label);
    } else {
        // set additional indel filters:
        const GatkVcfRecord& record(_recordBuffer[_indelIndex[0]]);
        const std::vector<std::string>& filters(record.GetFilter());
        const unsigned n_filt(filters.size());
        if((n_filt!=1) || filters[0] != "PASS") {
            rinfo.filters = filters;
        }
        
        // set additional minq:
        rinfo.qual.str=record.GetQual().c_str();
        rinfo.qual.is_valid=checked_double_parse(rinfo.qual.str,rinfo.qual.val);
        rinfo.gq.str=record.GetSampleVal("GQ");
        rinfo.gq.is_valid=checked_double_parse(rinfo.gq.str,rinfo.gq.val);

        _gti.clear();
        if(! record.GetGT().empty()) {
            parse_gt(record.GetGT().c_str(),_gti);
        }
        if(_gti.size() == 2) {
            if((_gti[0]==0 && _gti[1]>0) || (_gti[1]==0 && _gti[0]>0)){ rinfo.copyn=1; }
        }
    }


    // 2) modify site records according to overlapping filter status (or mark all as IndelConflict)
    //
    bool is_edit(true);
    std::vector<refedit> edits;
    for(unsigned record_index(0);record_index<n_records;++record_index) {
        GatkVcfRecord& record(_recordBuffer[record_index]);
        const int pos(record.GetPos());
        const bool is_in_indel((pos>=_bufferStartPos) && (pos<=_bufferEndPos));
        if(! is_in_indel) continue;
        const unsigned offset(pos-_bufferStartPos);
        adjust_overlap_record(_opt,rinfo,offset,record,is_edit,edits);
        // regroom record to account for quality value changes, etc:
        GroomInputRecord(record);
    }

    // 3) modify indel records according to any site conflicts or hemizygous snps present:

    // this is 90% done, but no easy way to make the per-allele tag adjustment reliable w/o parsing
    // header for all cases first. not worth pursuing for now... 
#if 0
    if(is_edit && (! edits.empty())) {
        // we should only get here for simple het deletions:
        GatkVcfRecord& record(_recordBuffer[_indelIndex[0]]);
        std::string allele(record.GetRef());
        bool is_diff(false);
        for(unsigned i(0);i<edits.size();++i) {
            if(allele[edits[i].first+1] != edits[i].second) {
                allele[edits[i].first+1] = edits[i].second;
                is_diff=true;
            }
        }
        if(is_diff) {
            // 1) insert new alternate allele
            // 2) update GT
            // 3) modify or delete all other allele dependent tags (this might just be AD in practice)
            std::vector<std::string>& alt(record.GetAlt());
            alt.insert(alt.begin(),allele);
            record.SetSampleVal("GT","1/2");
        }
    }
#endif


#ifdef VDEBUG 
    if(true) {
        std::cerr << "VDEBUG output: indel count: " << _indelIndex.size() << "\n";
        for(unsigned i(0);i<n_records;++i) {
            _recordBuffer[i].WriteUnaltered(std::cerr);
        }
    }
#endif
}



void
VcfRecordBlocker::
ProcessRecordBuffer() {

    if(_recordBuffer.empty()) return;

    // every record buffer should contain an indel:
    assert(_indelIndex.size() > 0);

    const unsigned n_records(_recordBuffer.size());

    // if there's only one record, assume this is a simple insertion and don't process
    // it for overlap information:
    if(n_records>1) GroomRecordBuffer();
    
    // send recordbuffer on for printing/blocking:
    for(unsigned i(0);i<n_records;++i) ProcessRecord(_recordBuffer[i]);
    _indelIndex.clear();
    _recordBuffer.clear();
}
