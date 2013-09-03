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

#include "blt_exception.hh"
#include "parse_util.hh"
#include "vcf_util.hh"

#include "boost/foreach.hpp"

#include <cassert>
#include <cctype>

#include <iostream>
#include <sstream>


struct gt_parse_helper {

    // return is_valid_genotype
    static
    bool
    start(const char* gt,
          std::vector<int>& gti,
          const bool is_badend) {
        gti.clear();
        if (isdigit(*gt)) return digit(gt,gti,is_badend);

        switch (*gt) {
        case '.' :  return unknown(gt,gti,is_badend);
        default: return false;
        }
    }

private:

    static
    bool
    unknown(const char* gt,
            std::vector<int>& gti,
            const bool is_badend) {
        gt++;
        gti.push_back(-1);
        switch (*gt) {
        case '\0' : return true;
        case '|' :
        case '/' : return sep(gt,gti,is_badend);
        default : return is_badend;
        }
    }

    static
    bool
    sep(const char* gt,
        std::vector<int>& gti,
        const bool is_badend) {
        gt++;
        if (isdigit(*gt)) return digit(gt,gti,is_badend);
        switch (*gt) {
        case '.' : return unknown(gt,gti,is_badend);
        default : return false;
        }
    }

    static
    bool
    digit(const char* gt,
          std::vector<int>& gti,
          const bool is_badend) {
        int val(0);
        while (isdigit(*gt)) {
            val = val*10 + static_cast<int>(*gt-'0');
            gt++;
        }
        gti.push_back(val);

        switch (*gt) {
        case '\0' : return true;
        case '|' :
        case '/' : return sep(gt,gti,is_badend);
        default : return is_badend;
        }
    }
};




void
parse_gt(const char* gt,
         std::vector<int>& gti,
         const bool is_allow_bad_end_char) {

    assert(NULL != gt);

    if (! gt_parse_helper::start(gt,gti,is_allow_bad_end_char)) {
        std::ostringstream oss;
        oss << "ERROR: can't parse genotype string: '" << gt << "'\n";
        throw blt_exception(oss.str().c_str());
    }
}



bool
is_variant_record(
    const char* const* word,
    std::vector<int>& gtparse) {

    const char* altstr(word[VCFID::ALT]);
    if (0==strcmp(".",altstr)) return false;

    parse_gt(get_format_string_nocopy(word,"GT"),gtparse,true);

    BOOST_FOREACH(const int allele, gtparse) {
        if (allele>0) return true;
    }
    return false;
}



void
get_vcf_end_record_range(
    const char* const* word,
    unsigned& begin_pos,
    unsigned& end_pos) {

    // get begin pos:
    const char* posstr(word[VCFID::POS]);
    begin_pos=(parse_unsigned(posstr));

    // get end pos:
    static const char* endkey = "END=";
    static const unsigned endsize = strlen(endkey);

    const char* endstr(strstr(word[VCFID::INFO],"END="));
    if (NULL==endstr) {
        end_pos = begin_pos;
    } else {
        endstr += endsize;
        end_pos = parse_unsigned(endstr);
    }
}



void
get_vcf_record_range(
    const char* const* word,
    unsigned& begin_pos,
    unsigned& end_pos) {

    get_vcf_end_record_range(word,begin_pos,end_pos);

    if (begin_pos != end_pos) return;

    // no END tag -- check to see if this is an indel record:
    const char* refStr(word[VCFID::REF]);
    const unsigned refLen(strlen(refStr));

    bool isIndel(false);
    {
        const char* altStr(word[VCFID::ALT]);

        // make sure there is an alternate:
        if (0!=strcmp(altStr,".")) {
            const char* tmp_ptr;
            while (NULL != (tmp_ptr=strchr(altStr,','))) {
                if ((tmp_ptr-altStr)!= refLen) isIndel=true;
                altStr = tmp_ptr+1;
            }
            if (strlen(altStr) != refLen) isIndel=true;
        }
    }
    if (! isIndel) return;

    // if it's an indel the first position isn't really "called"
    begin_pos += 1;
    end_pos += (refLen-1);

}
