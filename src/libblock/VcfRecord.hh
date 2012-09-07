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
#ifndef __VCF_RECORD_HH
#define __VCF_RECORD_HH

#include "istream_line_splitter.hh"
#include "string_util.hh"
#include "vcf_util.hh"

#include <cassert>
#include <cstring>

#include <iosfwd>



struct VcfRecord {

    VcfRecord(const istream_line_splitter& vparse);

    virtual ~VcfRecord() {}

    const std::string& GetChrom() const { return _chrom; }

    int GetPos() const { return _pos; }

    void SetPos(const int& pos) { _pos = pos; }

    const std::string& GetId() const { return _id; }

    const std::string& GetRef() const { return _ref; }

    std::vector<std::string>& GetAlt() { return _alt; }
    const std::vector<std::string>& GetAlt() const { return _alt; }

    const std::string& GetQual() const { return _qual; }

    void SetQual(const char* val) {
        if((NULL == val) || (*val=='\0')) {
            _qual = ".";
        } else {
            _qual = val;
        }
        IsSampleModified();
    }

    // this definition will pick up MNPs, we can refine this later if need be...
    bool IsIndel() const {
        if(GetRef().size() > 1) return true;
        const unsigned as(GetAlt().size());
        for(unsigned i(0);i<as;++i) {
            if(GetAlt()[i].size() > 1) return true;
        }
        return false;
    }

    // strictly, we need to check gt as well to determin if the site is actually genotyped as variant,
    //
    bool IsVariant() const {
        return (! GetAlt().empty());
    }

    bool IsNonvariantBlock() const {
        return ((GetRef().size() != 1) && (! IsVariant()));
    }

    const std::vector<std::string>& GetFilter() const { return _filt; }

    /// clear all filters and set to PASS state:
    void PassFilter() {
        _filt.clear();
        _filt.push_back("PASS");
    }

    /// add filter if it doesn't already exist:
    void AppendFilter(const char* val) {
        const unsigned fs(_filt.size());

        if((fs == 1) && (_filt[0] == "PASS")) {
            _filt.clear();
        }
        
        for(unsigned i(0);i<fs;++i) {
            if(_filt[i] == val) return;
        }
        _filt.push_back(val);
    }

    const char*
    GetInfoVal(const char* key) const {
        assert(NULL != key);
        const unsigned ic(_info.size());
        for(unsigned i(0);i<ic;++i){
            const size_t index(_info[i].find('='));
            if(index == std::string::npos) continue;
            if(0 == _info[i].compare(0,index,key)) {
                return _info[i].c_str()+index+1;
            }
        }
        return NULL;
    }

    void
    SetInfoVal(const char* key,
               const char* val) {
        assert(NULL != key);
        assert(NULL != val);
        const unsigned ic(_info.size());
        for(unsigned i(0);i<ic;++i){
            const size_t index(_info[i].find('='));
            if(index == std::string::npos) continue;
            if(0 == _info[i].compare(0,index,key)) {
                _info[i].replace(index,std::string::npos,val);
                return;
            }
        }
        _info.push_back(std::string(key)+"="+val);
    }

    void DeleteInfoKeyVal(const char* key) {
        const unsigned ic(_info.size());
        for(unsigned i(0);i<ic;++i){
            const size_t index(_info[i].find('='));
            if(index == std::string::npos) continue;
            if(0 == _info[i].compare(0,index,key)) {
                _info.erase(_info.begin()+i);
                return;
            }
        }
    }

    // client's responsibility to not insert repeats:
    void
    AppendInfo(const char* info) {
        _info.push_back(std::string(info));
    }

    void ClearInfo() { _info.clear(); }

    const char*
    GetSampleVal(const char* key) const {

        const unsigned fs(_format.size());
        for(unsigned i(0);i<fs;++i){
            if(_format[i] == key ) {
                return _sample[i].c_str();
            }
        }
        return NULL;
    }

    bool
    GetSampleValStr(const char* key,
                    std::string& val) const {

        const char* s(GetSampleVal(key));
        if(NULL==s) {
            val.clear();
            return false;
        } else {
            val = s;
            return true;
        }
    }

    void
    SetSampleVal(const char* key,
                 const char* val) {
 
        IsSampleModified();
        const unsigned fs(_format.size());
        for(unsigned i(0);i<fs;++i){
            if(_format[i] == key ) {
                _sample[i] = val;
                return;
            }
        }
        // add key if not found
        _format.push_back(key);
        _sample.push_back(val);
    }

    void
    DeleteSampleKeyVal(const char* key)
    {
        int deli = -1;
        const unsigned fs(_format.size());
        for(unsigned i(0);i<fs;++i){
            if(_format[i] == key ) {
                deli = i;
                break;
            }
        }

        if (deli < 0) return;
        _format.erase(_format.begin()+deli);
        _sample.erase(_sample.begin()+deli);
        IsSampleModified();
    }

    void 
    ReplaceSampleKey(const char* oldval,
                     const char* newval) {
        const unsigned fs(_format.size());
        for(unsigned i(0);i<fs;++i){
            if(_format[i] == oldval ) {
                _format[i] = newval;
                IsSampleModified();
                return;
            }
        }
    }

    void ClearSample() {
        _format.clear();
        _sample.clear();
        IsSampleModified();
    }

    void
    Write(const std::string& printChrom,
          const int printPos,
          const std::string& refPlaceholder,
          std::ostream& os) const;
        
    void
    Write(std::ostream& os) {
        static const std::string RPH = "..";
        const std::string* refvalptr = &(GetRef());

        // erase reference for non-variants across the board:
        if ((refvalptr->size() == 1) &&
            (! IsVariant())) {
            refvalptr = &RPH;
        }
        Write(GetChrom(), GetPos(), *refvalptr, os);
    }

    void
    WriteUnaltered(std::ostream& os) const {
        Write(GetChrom(), GetPos(), GetRef(), os);
    }

protected:

    virtual
    void IsSampleModified() {}

private:


    static
    void
    DumpVectorString(const std::vector<std::string>& v,
                     const char delimiter,
                     std::ostream& os);

    static void
    DumpAltString(const std::vector<std::string>& v,
                  std::ostream& os) {
        DumpVectorString(v,',',os);
    }

    static void
    DumpInfoString(const std::vector<std::string>& v,
                   std::ostream& os) {
        DumpVectorString(v,';',os);
    }

    static void
    DumpFormatString(const std::vector<std::string>& v,
                     std::ostream& os) {
        DumpVectorString(v,':',os);
    }

    static
    void 
    Splitter(const char* str,
             const char delimiter,
             std::vector<std::string>& v) {
        v.clear();
        if (NULL==str) return;
        if ('\0'==*str) return;
        if (0==strcmp(str,".")) return;
        split_string(str,delimiter,v);
    }

private:
    std::string _chrom;
    int _pos;
    std::string _id;
    std::string _ref;
    std::vector<std::string> _alt;
    std::string _qual;
    std::vector<std::string> _filt;
    std::vector<std::string> _info;
    std::vector<std::string> _format;
    std::vector<std::string> _sample;
};

//std::ostream& operator<<(std::ostream& os, const VcfRecord& vcfr);


#endif
