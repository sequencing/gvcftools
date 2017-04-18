// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2015 Illumina, Inc.
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
/// \author Chris Saunders and Subramanian Shankar Ajay
///


#include "related_sample_util.hh"
#include "string_util.hh"
#include "tabix_streamer.hh"
#include "vcf_util.hh"

#include "boost/algorithm/string/predicate.hpp"
#include "boost/foreach.hpp"

#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdlib>

#include <fstream>
#include <iostream>



namespace {
std::ostream& log_os(std::cerr);
}



std::ostream&
operator<<(std::ostream& os, const vcf_pos& vpos)
{
    os << "pos: " << vpos.pos << " is_indel: " << vpos.is_indel;
    return os;
}



static
bool
get_digt_code(const char* const* word,
              std::vector<int>& digt_code) {

    const char* gtstr(get_format_string_nocopy(word,"GT"));
    if (gtstr == NULL)
    {
        digt_code.clear();
        digt_code.push_back(-1);
    }
    else
    {
        parse_gt(gtstr,digt_code,true);
    }
    return (digt_code.size()==2 && digt_code[0]>=0 && digt_code[1]>=0);
}



bool
snp_type_info::
get_indel_allele(
    std::string& indel_ref,
    std::vector<std::string>& allele,
    const char* const* word) const
{
    allele.clear();

    indel_ref=(word[VCFID::REF]);
    const char* alt(word[VCFID::ALT]);

    std::vector<std::string> altWord;
    split_string(alt,',',altWord);

    get_digt_code(word,_gtcode);

    bool is_standard_diploid(true);

    BOOST_FOREACH(const int gt, _gtcode)
    {
        if (gt==0) {
            allele.push_back(indel_ref);
        }
        else if (gt<0)
        {
            allele.push_back("X");
        }
        else
        {
            assert(gt<=static_cast<int>(altWord.size()));
            allele.push_back(altWord[gt-1]);
        }
    }

    //return (is_standard_diploid && _gtcode.size()==2);
    return ((is_standard_diploid && _gtcode.size() == 2) || (_gtcode.size() == 1 && _gtcode[0] >= 0));
}



bool
snp_type_info::
get_site_allele(
    std::vector<char>& allele,
    const char* const* word,
    const unsigned offset,
    const char ref_base) const
{

    allele.clear();

    get_digt_code(word,_gtcode);

    bool is_standard_diploid(true);

    BOOST_FOREACH(const int gt, _gtcode)
    {
        if (gt==0) {
            allele.push_back(ref_base);
            continue;
        }

        const char* alt(word[VCFID::ALT]);
        for (int ai(0); (ai+1)<gt; alt++) {
            if (! *alt) break;
            if ((*alt)==',') ai++;
        }

        if ((! *alt) || (gt<0)) {
            allele.push_back('N');
            is_standard_diploid=false;
        } else {
            allele.push_back(alt[offset]);
        }
    }

    //return (is_standard_diploid && _gtcode.size()==2);
    return ((is_standard_diploid && _gtcode.size() == 2) || (_gtcode.size() == 1 && _gtcode[0] >= 0));
}


bool
snp_type_info::is_ref() const{

    BOOST_FOREACH(const int gt, _gtcode)
    {
        if(gt == 0)
            continue;
        else
            return false;
    }
    return true;
}                       


// extract unsigned value for key from the vcf info field
bool
snp_type_info::
get_info_unsigned(const char* info,
                  const char* key,
                  unsigned& val) {

    const char* tmp(NULL);
    do {
        if (NULL != tmp) info=tmp+1;
        if (0!=strncmp(info,key,strlen(key))) continue;
        if (NULL==(info=strchr(info,'='))) return false;
        const char* s(info+1);
        val=parse_unsigned(s);
        return true;
    } while (NULL != (tmp=strchr(info,';')));
    return false;
}



// extract float value for key from the vcf info field
bool
snp_type_info::
get_info_float(const char* info,
               const char* key,
               float& val) {

    const char* tmp(NULL);
    do {
        if (NULL != tmp) info=tmp+1;
        if (0!=strncmp(info,key,strlen(key))) continue;
        if (NULL==(info=strchr(info,'='))) return false;
        const char* s(info+1);
        val=parse_double(s);
        return true;
    } while (NULL != (tmp=strchr(info,';')));
    return false;
}



bool
snp_type_info::
get_format_float(const char* const* word,
                 const char* key,
                 float& val) {

    const char* str(get_format_string_nocopy(word,key));
    if (NULL==str) return false;
    if ('\0'==*str) return false;
    if ('.'==*str) {
        if ('\0'==*(str+1)) return false;
        if (':'==*(str+1)) return false;
    }
    val=parse_double(str);
    return true;
}


bool
snp_type_info::
get_format_unsigned(const char* const* word,
                    const char* key,
                    unsigned& val) {

    const char* str(get_format_string_nocopy(word,key));
    if (NULL==str) return false;
    val=parse_unsigned(str);
    return true;
}



site_crawler::
site_crawler(const sample_info& si,
             const unsigned /*sample_id*/,
             const shared_crawler_options& opt,
             const char* chr_region,
             const reference_contig_segment& ref_seg,
             const bool is_store_header,
             const bool is_return_indels)
    : _is_call(false)
    , _chrom(NULL)
    , _si(si)
    //, _sample_id(sample_id)
    , _opt(opt)
    , _chr_region(chr_region)
    , _is_return_indels(is_return_indels)
    , _tabs(NULL)
    , _is_sample_begin_state(true)
    , _is_sample_end_state(false)
    , _is_reference(false)
    , _next_file(0)
    , _ref_seg(ref_seg)
    , _n_word(0)
    , _locus_size(0)
    , _locus_offset(0)
    , _skip_call_begin_pos(0)
    , _skip_call_end_pos(0)
    , _is_site_allele_current(false)
    , _is_indel_allele_current(false)
{
    enum {MAX_WORD=50};
    memset(_prev_word,0,sizeof(char *)*(MAX_WORD));
    memset(_next_word,0,sizeof(char *)*(MAX_WORD));
    update(is_store_header);
}



site_crawler::
~site_crawler() {
    if (NULL != _tabs) delete _tabs;
    if (NULL != _prev_word[0]){
        unsigned idx(0);
        while (idx < _n_word)
            free(_prev_word[idx++]);
    }
}




void
site_crawler::
dump_state(std::ostream& os) const {
    const std::string& afile(_si.file);
    os << "LOCUS_CRAWLER STATE:\n";
    os << "\tchrom: " << chrom();
    os << "\tposition: " << pos() << " offset: " << _locus_offset << "\n";
    os << "\t_is_call: " << _is_call << "\n";
    os << "\t_is_site_allele_current: " << _is_site_allele_current;
    os << "\t_is_sample_begin_state: " << _is_sample_begin_state;
    os << "\t_is_sample_end_state: " << _is_sample_end_state;
    os << "\tis_indel: " << is_indel() << "\n";
    os << "\tfile: '" << afile << "'\n";
    os << "\tline: '";
    dump_line(os);
    os << "'\n";
}





bool
site_crawler::
update_site_allele() const
{
    const char ref_base(_ref_seg.get_base(pos()-1));
    const bool retval(_opt.sti().get_site_allele(_site_allele,_word,_locus_offset,ref_base));
    _is_site_allele_current=true;
    return retval;
}



bool
site_crawler::
update_indel_allele() const
{
    const bool retval(_opt.sti().get_indel_allele(_indel_ref,_indel_allele, _word));
    _is_indel_allele_current=true;
    return retval;
}



bool
site_crawler::
update_allele() const {
    if (is_indel()) {
        return update_indel_allele();
    } else {
        return update_site_allele();
    }
}



static const char sep('\t');


// return true if current position in record is valid and usable
bool
site_crawler::
process_record_line(char* line)
{
    static const unsigned MAX_WORD(50);

    // do a low-level tab parse:
    {
        char* p(line);
        _word[0]=p;
        _n_word=1;
        while (true) {
            if ((*p == '\n') || (*p == '\0')) break;
            if (*p == sep) {
                *p = '\0';
                _word[_n_word++] = p+1;
                if (_n_word == MAX_WORD) break;
            }
            ++p;
        }
        // allow for optional extra columns in each file format:
        if (_n_word<_opt.sti().col_count()) {
            log_os << "ERROR: Consensus record has " << _n_word << " column(s) but expecting at least " << _opt.sti().col_count() << "\n";
            dump_state(log_os);
            exit(EXIT_FAILURE);
        }
    }

    const vcf_pos last_vpos(vpos());
    _chrom=_opt.sti().chrom(_word);
    _vpos.pos=_opt.sti().pos(_word);

    _is_site_allele_current=false;
    _is_indel_allele_current=false;

    _vpos.is_indel=(_opt.sti().get_is_indel(_word));

    if (pos()<1) {
        log_os << "ERROR: gvcf record position less than 1. position: " << pos() << " ";
        dump_state(log_os);
        exit(EXIT_FAILURE);
    }

    if (_opt.is_region()) {
        // deal with vcf records after the region of interest:
        if (pos()>_opt.region_end) {
            _is_sample_begin_state = false;
            _is_sample_end_state = true;
            return true;
        }
    } else {
        if (pos()>static_cast<pos_t>(_ref_seg.end())) {
            log_os << "ERROR: allele file position exceeds final position in reference sequence segment . position: "
                   << pos() << " ref_contig_end: " << _ref_seg.end() << "\n";
            dump_state(log_os);
            exit(EXIT_FAILURE);
        }
    }

    if (! _opt.sti().get_nonindel_ref_length(pos(),is_indel(),_word,_locus_size)) {
        //log_os << "ERROR: failed to parse locus at pos: "  << pos << "\n";
        log_os << "WARNING: failed to parse locus at: "  << vpos() << "\n";
        dump_state(log_os);
        //exit(EXIT_FAILURE);
        _locus_size=0;
    }

    _locus_offset=0;

    // deal with vcf records which fully proceed the region of interest:
    if (_opt.is_region()) {
        if ((pos()+_locus_size-1)<_opt.region_begin) return false;
    }

    //const bool last_is_call(is_call);
    _is_call = _opt.sti().get_is_call(_word,pos(),_skip_call_begin_pos,_skip_call_end_pos);

    _n_total = _opt.sti().total(_word);

    if (is_indel()) {
        if (! _is_return_indels)
        {
            _vpos.pos=last_vpos.pos;
            _locus_size=0;
            return false;
        }
        else
        {
            _locus_size=0;
        }
    }

    if (is_any_call()) {
        _is_call=update_allele();
    }

    // don't allow failed block read-through, so that we can get through indel-overlap errors
    //if(! is_call) {
    //    if(_locus_size>1) _locus_size=1;
    //}

    if (! _is_sample_begin_state) {
        //if (! (last_vpos < vpos()) ) {
        if (vpos() < last_vpos ) {
            if (_opt.is_murdock_mode) {
                _vpos=last_vpos;
                _locus_size=0;
                return false;
            } else {
                log_os << "ERROR: unexpected position order in variant file. current_pos: "
                       << pos() << " last " << last_vpos << "\n";
                dump_state(log_os);
                exit(EXIT_FAILURE);
            }
        }
    } else {
        _is_sample_begin_state=false;
    }

    // deal with vcf records which partially overlap the region of interest:
    if (_opt.is_region()) {
        if (pos()<_opt.region_begin) return false;
    }

    return true;
}



void
site_crawler::
dump_header(std::ostream& os) const {

    const unsigned header_size(_header.size());
    for (unsigned i(0); (i+1)<header_size; ++i) {
        os << _header[i] << '\n';
    }
}


bool
site_crawler::
rewind_site(const vcf_pos& lpos) {

    if(_locus_offset > 0){
        while(_vpos.pos > lpos.pos){
            _locus_offset--;
            _vpos.pos--;
        }
        _is_call = update_allele();
        _is_reference = _opt.sti().is_ref();
    }
    else{
        unsigned idx(0);
        while(idx < _n_word){
            _next_word[idx] = _word[idx];
            _word[idx] = _prev_word[idx];
            idx++;
        }
        _vpos.pos = _opt.sti().pos(_word);
        _is_site_allele_current = false;
        _is_indel_allele_current = false;
        _vpos.is_indel = (_opt.sti().get_is_indel(_word));
        _opt.sti().get_nonindel_ref_length(pos(),is_indel(),_word,_locus_size);

        if(is_indel())
            _locus_size = 1;
            
        _locus_offset = _locus_size - 1;
        _vpos.pos += _locus_offset;
        while(_vpos.pos > lpos.pos){
            _locus_offset--;
            _vpos.pos--;
        }
        _is_call = update_allele();
        _is_reference = _opt.sti().is_ref();
    }

    return true;
}


void
site_crawler::
update(bool is_store_header) {
    if (_is_sample_end_state) return;


    // move on to a new locus:
    while (true) {
        // continue crawling through a multibase locus:
        //
        // at present (201202) this only means compressed block
        // entries in gVCF -- long indels should be skipped over
        // rather than walked
        //
        _locus_offset++;
        if (_locus_offset < _locus_size) {
            // check for pos moving past the end of region of interest:
            if (_opt.is_region()) {
                if ((pos()+1)>(_opt.region_end)) {
                    _is_sample_begin_state = false;
                    _is_sample_end_state = true;
                    return;
                }
            }

            _vpos.pos++;
            _is_site_allele_current = false;

            if (_opt.is_region()) {
                // check for pos preceding the start of region of interest in a multi-base record:
                if (pos()<_opt.region_begin) {
                    continue;
                }
            }

            if (pos()>static_cast<pos_t>(_ref_seg.end())) {
                log_os << "ERROR: allele file position exceeds final position in reference sequence segment. position: "
                       << pos() << " ref_contig_end: " << _ref_seg.end() << "\n";
                dump_state(log_os);
                exit(EXIT_FAILURE);
            }

            if (is_site_call()) {
                const char ref_base=_ref_seg.get_base(pos()-1);
                if (! _opt.sti().get_site_allele(_site_allele,_word,_locus_offset,ref_base)) {
                    log_os << "ERROR: Failed to read site genotype from record:\n";
                    dump_state(log_os);
                    exit(EXIT_FAILURE);
                }
            }
            return;
        }

        // start new/next file:
        if (NULL == _tabs) {
            if (_next_file >= 1) {
                _is_sample_begin_state = false;
                _is_sample_end_state = true;
                return;
            }
            const std::string& afile(_si.file);
            if (0 == _next_file) {
                // open a separate header streamer to get sample_name and optional header capture:
                tabix_header_streamer ths(afile.c_str());
                while (ths.next()) {
                    const char* line(ths.getline());
                    if (NULL == line) break;
                    if (is_store_header) {
                        _header.push_back(line);
                    }
                    if (_sample_name.empty() && boost::starts_with(line,"#CHROM")) {
                        std::vector<std::string> words;
                        split_string(line,'\t',words);
                        if (words.size()>VCFID::SAMPLE) {
                            _sample_name = words[VCFID::SAMPLE];
                        } else {
                            _sample_name = "UNKNOWN";
                        }
                    }
                }
            }
            _tabs=new tabix_streamer(afile.c_str(),_chr_region);
            if (! _tabs) {
                log_os << "ERROR:: Can't open gvcf file: " << afile << "\n";
                exit(EXIT_FAILURE);
            }
            _next_file++;
        }

        if (_next_word[0]){
            unsigned idx(0);
            while(idx < _n_word){
                _word[idx] = _next_word[idx];
                _next_word[idx] = NULL;
                idx++;
            }
            _vpos.pos = _opt.sti().pos(_word);
            _locus_offset = 0;
            _is_site_allele_current = false;
            _is_indel_allele_current = false;
            _vpos.is_indel=(_opt.sti().get_is_indel(_word));
            _opt.sti().get_nonindel_ref_length(pos(),is_indel(),_word,_locus_size);
            _is_call=update_allele();

            return;
        }

        unsigned idx(0);
        if ( ! _is_sample_begin_state){
            while (idx < _n_word){
                size_t size = strlen(_word[idx]) + sizeof(char) ;

                if (NULL != _prev_word[idx])
                    free(_prev_word[idx]);
                _prev_word[idx] = (char *)malloc(size);
                strncpy(_prev_word[idx], _word[idx], size);
                idx++;
            }
        }


        // read through file to get to a data line:
        bool is_eof(true);
        char* line(NULL);
        while (_tabs->next()) {
            line=_tabs->getline();
            if (NULL == line) break;
            assert(strlen(line));
            if (line[0] == '#') continue;
            is_eof=false;
            break;
        }

        if (! is_eof) {
            const bool is_valid=process_record_line(line);
            if (is_valid) return;
        } else {
            // if eof, go to next file or terminate:
            delete _tabs;
            _tabs=NULL;
        }
    }
}



void
site_crawler::
dump_line(std::ostream& os) const {
    if (_is_sample_end_state) return;
    for (unsigned i(0); i<_n_word; ++i) {
        if (i) os << sep;
        os << _word[i];
    }
}



pos_reporter::
pos_reporter(const std::string& filename,
             const std::vector<std::string>& sample_label)
    : is_pos_report(false),
      pos_fs_ptr(0),
      sample_size(sample_label.size()),
      _sample_label(sample_label) {

    if (filename.empty()) return;
    pos_fs_ptr = new std::ofstream(filename.c_str());
    if (! pos_fs_ptr) {
        log_os << "ERROR: Failed to allocate stream pointer for file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
    if (! *pos_fs_ptr) {
        log_os << "ERROR: Failed to open output file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
    is_pos_report=true;
}



pos_reporter::
~pos_reporter() {
    if (is_pos_report) delete pos_fs_ptr;
}



void
pos_reporter::
print_pos(const site_crawler* sa) {
    if (! is_pos_report) return;
    *pos_fs_ptr << "EVENT\t" << sa[0].chrom() << "\t" << sa[0].pos() << "\n";
    for (unsigned i(0); i<sample_size; ++i) {
        *pos_fs_ptr << _sample_label[i] << "\t";
        *pos_fs_ptr << "\t";
        sa[i].dump_line(*pos_fs_ptr);
        *pos_fs_ptr << "\n";
    }
}

