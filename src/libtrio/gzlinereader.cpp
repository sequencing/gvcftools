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

#include "gzlinereader.hh"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <iostream>



namespace {
std::ostream& log_os(std::cerr);
}



gzlinereader::
gzlinereader(const char* filename)
  : _filename(filename),
    _lineno(1),
    _buf(new char[_bufsize+1]),
    _start(0),
    _end(0),
    _is_last(false) {

    _zfp=gzopen(filename,"rb");
    if(_zfp==Z_NULL) {
        log_os << "ERROR: failed to open file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }
}



gzlinereader::
~gzlinereader() {
    gzclose(_zfp);
    delete _buf;
}



char*
gzlinereader::
getline() {
    unsigned p(_start);
    while(true) {
        while(p==_end) {
            if(_is_last) return NULL;
            read_buffer();
            p=_start;
        }

        for(;p<_end;++p){
            if(_buf[p]!='\n') continue;
            _lineno++;
            char* retptr(_buf+_start);
            _start=p+1;
            _buf[p]='\0'; 
            return retptr;
        }

        if(_is_last) {
            log_os << "ERROR: file is not terminated by newline: " << _filename << " line_no: " << _lineno << "\n";
            exit(EXIT_FAILURE);
        }

        if(_start==0) {
            log_os << "ERROR: line exceeds buffer size in file: " << _filename << " line_no: " << _lineno << "\n";
            exit(EXIT_FAILURE);
        }
    }
}



void
gzlinereader::
read_buffer() {
#ifdef GLR_DEBUG
   log_os << "\n...READING... start: " << _start << " end: " << _end << " is_last: " << _is_last << "\n"
          << "buf: XXX" << std::string(_buf,_end) << "XXX\n";
#endif
   assert(_start <= _end);
   assert((_start != 0) or (_end == 0));

   if(_is_last) {
       _start=0;
       _end=0;
       return;
   }

   const unsigned offset(_end-_start);
   if(offset != 0) {
       memmove(_buf,_buf+_start,offset);
   }

   const unsigned target_ret(_bufsize-offset);
      
   const int ret(gzread(_zfp,_buf+offset,target_ret));
   if(ret < 0) {
       int errnum;
       log_os << "ERROR: gzip error: " << gzerror(_zfp,&errnum) << "\n";
       exit(EXIT_FAILURE);
   }

   if(static_cast<unsigned>(ret)!=target_ret) _is_last=true; 
   _start=0;
   _end=ret+offset;

#ifdef GLR_DEBUG
   log_os << "...END_READING... start: " << _start << " end: " << _end << " is_last: " << _is_last << "\n"
          << "buf: XXX" << std::string(_buf,_end) << "XXX\n";
#endif
}


void
test() {
    gzlinereader glr("test.gz");

    const char* s;
    while(NULL != (s = glr.getline())) { log_os << s << "\n"; }
 
#if 0
   gzFile gz;
   gz = gzopen("test", "rb");
   assert( gz != NULL );

   static const unsigned bufsize(64);
   char* buf(new char[bufsize]);

   while(true) {
       int ret=gzread(gz,buf,bufsize);
       if(ret < 0) {
         int errnum;
         std::cerr << "ERROR: gzip error: " << gzerror(gz,&errnum) << "\n";
         exit(EXIT_FAILURE);
       }
       for(unsigned i(0);i<ret;++i){
         std::cerr << buf[i];
       }
       if(ret!=bufsize) break;
   }

  delete buf;
#endif
}

