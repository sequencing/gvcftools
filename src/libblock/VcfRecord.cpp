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
#include "VcfRecord.hh"

#include <iostream>
#include <sstream>


VcfRecord::
VcfRecord(const istream_line_splitter& vparse)
{
    const unsigned ws(vparse.n_word());
    if (static_cast<int>(ws) <= VCFID::INFO) {
        std::ostringstream oss;
        oss << "Too few fields (" << ws << ") in vcf record input:\n";
        vparse.dump(oss);
        throw blt_exception(oss.str().c_str());
    }

    _chrom = vparse.word[VCFID::CHROM];

    const char* pos_ptr(vparse.word[VCFID::POS]);
    _pos = parse_unsigned(pos_ptr);

    _id = vparse.word[VCFID::ID];

    _ref = vparse.word[VCFID::REF];
    assert(_ref.size() > 0);

    Splitter(vparse.word[VCFID::ALT],',',_alt);

    for (unsigned i(0); i<_alt.size(); ++i) {
        assert(_alt[i].size() > 0);
    }

    _qual = vparse.word[VCFID::QUAL];

    Splitter(vparse.word[VCFID::FILT],';',_filt);

    Splitter(vparse.word[VCFID::INFO],';',_info);

    if (ws > VCFID::FORMAT) {
        Splitter(vparse.word[VCFID::FORMAT],':',_format);
    }

    if (ws > VCFID::SAMPLE) {
        Splitter(vparse.word[VCFID::SAMPLE],':',_sample);
    }

    // by the vcf spec, we can drop trailing fields in any sample:
    if (_sample.size() < _format.size())
    {
        _sample.resize(_format.size(),".");
    }

    if (_format.size() != _sample.size()) {
        std::ostringstream oss;
        oss << "FORMAT and SAMPLE fields do not agree for vcf record:\n";
        vparse.dump(oss);
        throw blt_exception(oss.str().c_str());
    }


}


void
VcfRecord::
Write(const std::string& printChrom,
      const int printPos,
      const std::string& refPlaceholder,
      std::ostream& os) const {

    os << printChrom << '\t'
       << printPos << '\t'
       << _id << '\t'
       << refPlaceholder << '\t';

    DumpAltString(_alt,os);
    os << '\t'
       << _qual << '\t';

    DumpInfoString(_filt,os);
    os << '\t';
    DumpInfoString(_info,os);
    os << '\t';
    DumpFormatString(_format,os);
    os << '\t';
    DumpFormatString(_sample,os);
    os << '\n';
}



void
VcfRecord::
DumpVectorString(const std::vector<std::string>& v,
                 const char delimiter,
                 std::ostream& os) {
    if (v.empty()) {
        os << '.';
        return;
    }
    const unsigned vs(v.size());
    for (unsigned i(0); i<vs; ++i) {
        if (i) os << delimiter;
        os << v[i];
    }
}

