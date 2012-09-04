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
/// an efficient (and slightly unsafe) class for basic tab-delimited files, etc...
///

/// \author Chris Saunders
///
#ifndef ISTREAM_LINE_SPLITTER_HH__
#define ISTREAM_LINE_SPLITTER_HH__

#include <iosfwd>


struct istream_line_splitter {

    istream_line_splitter(std::istream& is,
                          const unsigned line_buf_size=8*1024,
                          const char word_seperator='\t',
                          const unsigned max_word=0)
        : _is(is)
        , _line_no(0)
        , _n_word(0)
        , _buf_size(line_buf_size)
        , _sep(word_seperator)
        , _max_word(max_word)
        , _buf(new char[_buf_size]) {

        if((0==_max_word) || (MAX_WORD_COUNT < _max_word)){
            _max_word=MAX_WORD_COUNT;
        }
    }

    ~istream_line_splitter() { if(NULL!=_buf) { delete [] _buf; _buf=NULL;} }

    unsigned
    n_word() const { return _n_word; }

    void
    dump(std::ostream& os) const;

    /// returns false for regular end of input:
    bool
    parse_line();


    enum { MAX_WORD_COUNT = 50 };
    char* word[MAX_WORD_COUNT];
private:
    std::istream& _is;
    unsigned _line_no;
    unsigned _n_word;
    unsigned _buf_size;
    char _sep;
    unsigned _max_word;
    char* _buf;
};



#if 0
{ //usage example:
    istream_line_splitter dparse(data_is);

    while(dparse.parse_line()) {
        static const unsigned col_count(46);
        if(dparse.n_word()!=col_count){
            std::ostringstream oss;
            oss << "ERROR: unexpected number of columns in paired export line:\n\n";
            dparse.dump(oss);
            throw blt_exception(oss.str().c_str());
        }
        
        for(unsigned i(1);(i+1)<col_count;++i){
            dparse.word[i][strlen(dparse.word[i])] = sep;
        }
        const char* nocompress_segment(dparse.word[0]);
        const char* compress_segment(dparse.word[1]);

        /// ....etc
    }
}
#endif


#endif
