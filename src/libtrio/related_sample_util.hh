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
/// \author Chris Saunders
///

#ifndef __RELATED_SAMPLE_UTIL_HH
#define __RELATED_SAMPLE_UTIL_HH

#include "parse_util.hh"
#include "pos_type.hh"
#include "reference_contig_segment.hh"

#include "boost/lexical_cast.hpp"

#include <cassert>
#include <cstdio>
#include <cstring>

#include <iosfwd>
#include <string>
#include <vector>



inline
double
ratio(const unsigned n, const unsigned d){
    return ((d==0)? 0. : (static_cast<double>(n)/static_cast<double>(d)));
}



struct site_stats_core {

    site_stats_core()
        : ref_size(0)
        , known_size(0)
        , some_mapped(0)
        , some_called(0)
        , all_mapped(0)
        , all_called(0)
        , incorrect(0)
        , snp_correct(0)
    {}

    unsigned ref_size;
    unsigned known_size;
    unsigned some_mapped;
    unsigned some_called;
    unsigned all_mapped;
    unsigned all_called;
    unsigned incorrect;
    unsigned snp_correct;
};



struct snp_param {

    snp_param()
        : min_gqx(0)
        , min_qd(0)
        , min_pos_rank_sum(0)
        , is_min_qd(false)
        , is_min_pos_rank_sum(false)
    {}

    double min_gqx;
    double min_qd;
    double min_pos_rank_sum;

    bool is_min_qd;
    bool is_min_pos_rank_sum;
};



struct snp_type_info {

    snp_type_info(const snp_param& sp,
                  const unsigned ccol,
                  const unsigned pcol) 
        : _sp(sp)
        , _chromcol(ccol)
        , _poscol(pcol)
    {}

    virtual
    ~snp_type_info() {}

    virtual
    bool
    get_is_call(char** word,
                const pos_t pos,
                const bool is_indel,
                pos_t& skip_call_begin_pos,
                pos_t& skip_call_end_pos) const = 0;

    virtual
    bool
    get_allele(std::pair<char,char>& allele,
               const char * const * word,
               const unsigned offset,
               const char ref_base) const = 0;

    virtual
    const char*
    score(const char * const * word) const = 0;

    virtual
    unsigned
    total(const char * const * word) const = 0;

    virtual
    unsigned
    col_count() const = 0; // total number of columns expected in a data line

    const char*
    chrom(const char* const * word) const {
        return word[_chromcol];
    }

    pos_t
    pos(const char* const * word) const {
        return boost::lexical_cast<pos_t>(word[_poscol]);
    }


    // Size of the site or multisite locus in reference bases. This has been 
    // assumed to be 1 until vcf support was added. On indel records this
    // returns 0.
    virtual
    bool
    get_nonindel_ref_length(const pos_t /*pos**/, const bool /*is_indel*/, const char * const * /*word*/,unsigned& result) const {
        result=1;
        return true;
    }
   
    virtual
    bool
    get_is_indel(const char * const * word) const = 0; 

protected:
    const snp_param& _sp;

private:
    const unsigned _chromcol;
    const unsigned _poscol;
};



struct snp_type_info_vcf : public snp_type_info {

    snp_type_info_vcf(const snp_param& sp) 
        : snp_type_info(sp,0,1)
    {}

    bool
    get_is_call(char** word,
                const pos_t pos,
                const bool is_indel,
		pos_t& skip_call_begin_pos,
		pos_t& skip_call_end_pos) const {

        update_skip_call_range(word,pos,is_indel,skip_call_begin_pos,skip_call_end_pos);
        bool is_bad_call(is_indel || (0!=strcmp(word[FILT_COL],"PASS")));
        if(is_bad_call) return false;
        if(_sp.min_gqx > 0) {
            float gqx(0);
            if(get_format_float(word,"GQX",gqx)) {
                if(gqx<_sp.min_gqx) return false;
            }
        }
#if 0
        if(_sp.vcf_min_ref_qual > 0) {
            const char* qword(word[QVAR_COL]);
            if(strcmp(".",qword)!=0) {
                unsigned digt_code[2];
                if(get_digt_code(word, digt_code) &&
                   (digt_code[0]==0) &&
                   (digt_code[1]==0)) {
                    float qual(casava::blt_util::parse_double(qword));
                    if(qual<_sp.vcf_min_ref_qual) return false;
                }
            }
        }
#endif
        if(_sp.is_min_qd) {
            float qd(0);
            if(get_info_float(word[INFO_COL],"QD",qd)) {
                if(qd<_sp.min_qd) return false;
            }
        }
        if(_sp.is_min_pos_rank_sum) {
            float pos_rank_sum(0);
            if(get_info_float(word[INFO_COL],"BaseQRankSum",pos_rank_sum)) {
                if(pos_rank_sum<_sp.min_pos_rank_sum) return false;
            }
        }
        return (pos<skip_call_begin_pos || pos>=skip_call_end_pos);
    }
  

    bool
    get_allele(std::pair<char,char>& allele,
               const char * const * word,
               const unsigned offset,
               const char ref_base) const;

    const char*
    score(const char * const * /*word*/) const {
        return "0";
    }

    unsigned
    total(const char * const * word) const {
        float dp;
        if(! get_format_float(word,"DP",dp)) return 0;
        return static_cast<unsigned>(dp);
    }

    unsigned
    col_count() const { // total number of columns expected in a data line
        return 10;
    }

    // Size of the site or multisite locus in reference bases. This has been 
    // assumed to be 1 until vcf support was added. On indel records this
    // returns 0.

    // returns false on error
    bool
    get_nonindel_ref_length(const pos_t pos, const bool is_indel, const char * const * word, unsigned& result) const { 
        result=0;
        if(is_indel) return true;
        
        const unsigned reflen(strlen(word[REF_COL]));
        unsigned iend;
        if(! get_info_unsigned(word[INFO_COL],"END",iend)) {
            result=reflen;
        } else {
            if( ! (iend>=pos && (reflen <= ((iend+1)-pos))) ) return false; 
            result = (iend+1)-pos;
        }
        return true;
    }

    // this detects both indels and unequal or equal length 'block-substutitons'
    bool
    get_is_indel(const char * const * word) const {
        const char* alt(word[ALT_COL]);
        const char* tmp_ptr;
        // no alternate:
        if(0==strcmp(alt,".")) return false;
        // breakend:
        if(NULL != (tmp_ptr=strchr(alt,'.'))) return true;
        const unsigned reflen(strlen(word[REF_COL]));
        // if alt is not '.' and reflen > 1, then this must be some sort of indel/subst:
        // pathological case is alt=".,." ... don't worry about that one.
        if(reflen>1) return true;
        // after the above test, we're essentially just looking for an alt with
        // length greater than one, indicating an insertion:
        while(NULL != (tmp_ptr=strchr(alt,','))){
            if((tmp_ptr-alt)!=static_cast<pos_t>(reflen)) return true;
            alt = tmp_ptr+1;
        } 
        return (strlen(alt)!=reflen);
    }

private:
    static
    char*
    get_format_string(const char* const * word, 
                      const char* key);

    static
    bool
    get_digt_code(const char * const * word,
                  unsigned digt_code[2]);


    static
    bool
    get_info_float(const char* info,
                   const char* key,
                   float& val);

    static
    bool
    get_info_unsigned(const char* info,
                      const char* key,
                      unsigned& val);

    static
    bool
    get_format_float(const char* const * word, 
                     const char* key,
                     float& val);

    // Size of the locus in reference bases. This has been assumed to
    // be 1 until vcf support was added.
    unsigned
    get_ref_length(const char * const * word) const { 
        return strlen(word[REF_COL]);
    }

    
    // if the current locus is an indel, mark out its range as no-call:
    //
    void
    update_skip_call_range(const char * const * word,
                           const pos_t pos,
                           const bool is_indel,
			   pos_t& skip_call_begin_pos,
			   pos_t& skip_call_end_pos) const {
        const unsigned locus_size=get_ref_length(word);
        if(locus_size<=1) return;
	if(! is_indel) return;
        // follow CASAVA begin-end convention: 0-indexed and end goes 1 past range
        const pos_t indel_begin_pos(pos+1);
        const pos_t indel_end_pos(pos+locus_size);
        if(skip_call_end_pos<pos){
            skip_call_begin_pos=indel_begin_pos;
        }
        skip_call_end_pos=std::max(skip_call_end_pos,indel_end_pos);
    }

    enum {
        POS_COL = 1,
        REF_COL = 3,
        ALT_COL = 4,
        QVAR_COL = 5,
        FILT_COL = 6,
        INFO_COL = 7,
        FORMAT_COL = 8,
        SAMPLE_COL = 9,
    };
};




// used to be a big stuct!!
struct sample_info {
    std::string file;
};


struct shared_crawler_options {
    shared_crawler_options()
        : stip(NULL)
        , is_murdock_mode(false)
    {}

    const snp_type_info& sti() const { return *stip; }

    bool is_region() const { return (! region.empty()); }

    snp_type_info* stip;
    std::string region;
    int region_begin;
    int region_end;
    bool is_murdock_mode;
};



struct tabix_streamer;


struct site_crawler {
    
    site_crawler(const sample_info& si,
                 const unsigned sample_id,
                 const shared_crawler_options& opt,
                 const char* chr_region,
                 const reference_contig_segment& ref_seg);

    ~site_crawler();

    void
    update();

    bool
    is_pos_valid() const { return (! _is_sample_end_state); }

    void
    dump_line(std::ostream& os) const;

    void
    dump_state(std::ostream& os) const;

    const char*
    chrom() const {
        static const char* unknown = "unknown";
        if(NULL == _chrom) return unknown;
        return _chrom;
    }

    enum {
        buff_size = 5000
    };

    pos_t pos;
    bool is_call;
    std::pair<char,char> allele;
    char score[buff_size];
    unsigned n_total;

private:

    bool
    process_record_line(char* line);

    const char* _chrom;
    const sample_info& _si;
    const unsigned _sample_id; // only useed for debuging...
    const shared_crawler_options& _opt;
    const char* _chr_region;

    tabix_streamer* _tabs;
    bool _is_sample_begin_state;
    bool _is_sample_end_state;
    unsigned _next_file;
    const reference_contig_segment& _ref_seg;

    enum { MAX_WORD=50 };
    char* _word[MAX_WORD];
    unsigned _n_word;

    unsigned _locus_size;
    unsigned _locus_offset;

    mutable pos_t _skip_call_begin_pos;
    mutable pos_t _skip_call_end_pos;
};



struct pos_reporter {

    pos_reporter(const std::string& filename,
                 const std::vector<std::string>& sample_label);

    ~pos_reporter();

    void
    print_pos(const site_crawler* sa);

private:

    bool is_pos_report;
    std::ofstream* pos_fs_ptr;
    unsigned sample_size;
    const std::vector<std::string>& _sample_label;
};

#endif
