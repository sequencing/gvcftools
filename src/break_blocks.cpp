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

#include "blt_exception.hh"
#include "compat_util.hh"
#include "gvcftools.hh"
#include "parse_util.hh"
#include "ref_util.hh"
#include "region_util.hh"
#include "stringer.hh"
#include "VcfHeaderHandler.hh"
#include "VcfRecord.hh"


#include "boost/program_options.hpp"

//#include <ctime>

#include <iostream>
#include <string>


namespace {
std::ostream& log_os(std::cerr);
}

std::string cmdline;


struct BreakOptions {

    BreakOptions()
        : outfp(std::cout)
    {}

    std::ostream& outfp;
    std::string ref_seq_file;
    region_util::region_t regions;
};



// process each vcf record for haploid setting:
//
struct BreakVcfRecordHandler {

    BreakVcfRecordHandler(const BreakOptions& opt)
        : _opt(opt)
        , _scp(opt.ref_seq_file.c_str())
    {}

    void
    process_line(const istream_line_splitter& vparse) {
        const unsigned nw(vparse.n_word());

        if(nw+1 != VCFID::SAMPLE) {
            log_os << "ERROR: unexpected number of fields in vcf record:\n";
            vparse.dump(log_os);
            exit(EXIT_FAILURE);
        }

        // 1. check if record region is in a target region at all,
        // if true an iterator provides begin,end ranges on
        // successive calls for the record region.
        //
        if(! is_record_in_region(vparse)) {
            vparse.write_line(_opt.outfp);
            return;
        }

        VcfRecord vcfr(vparse);

        bool is_deblock;
        unsigned end;
        while(get_next_record_region_interval(is_deblock,end)){
            VcfRecord vcfr2(vcfr);
            process_block(vcfr2,is_deblock,end);
            vcfr.SetPos(end+1);
            vcfr.SetRef(_scp.get_char(vcfr.GetChrom().c_str(),static_cast<int>(end+1)));
        }
        process_block(vcfr,is_deblock,end);
   }

private:
    bool
    is_record_in_region(const istream_line_splitter& vparse) {
        // determine if chromosome is new:
        if(_last_chrom.empty() || (0 != strcmp(_last_chrom.c_str(),vparse.word[0]))) {
            _last_chrom=vparse.word[0];
            const region_util::region_t::const_iterator i(_opt.regions.find(_last_chrom));
            _is_skip_chrom=((i==_opt.regions.end()) || (i->second.empty()));

            // setup region iterators:
            if(! _is_skip_chrom) {
                _rhead=i->second.begin();
                _rend=i->second.end();
            }
        }

        if(! _is_skip_chrom) {
            // get start pos:
            const char* posstr(vparse.word[1]);
            _begin_pos=(parse_unsigned(posstr));

            // get end pos:
            static const char* endkey = "END=";
            static const unsigned endsize = strlen(endkey);
            const char* endstr(strstr(vparse.word[VCFID::INFO],"END="));
            if(NULL==endstr) {
                _end_pos = _begin_pos;
            } else {
                endstr+=endsize;
                _end_pos=parse_unsigned(endstr);
            }

            while(_rhead != _rend) {
                if(_begin_pos>_rhead->second) {
                    _rhead++;
                } else {
                    return(_end_pos>_rhead->first);
                }
            }
            _is_skip_chrom=true;
        }
        return false;
    }

    // if is_record_in_region is true, call this function to get intercepting intervals
    // returns false when no more intervals exist
    bool
    get_next_record_region_interval(bool& is_in_region,
                                    unsigned& end) {

        assert(_begin_pos <= _end_pos);

        if(_begin_pos > _rhead->second) _rhead++;
        if(_rhead == _rend) { // no haploid regions left
            is_in_region = false;
            end = _end_pos;
            return false;
        }

        // our next interval:
        if(_begin_pos <= _rhead->first) {
            end = std::min( _end_pos,_rhead->first);
        } else {
            end = std::min( _end_pos, _rhead->second);
        }

        // test for intercept:
        is_in_region = ((_begin_pos <= _rhead->second) && (end > _rhead->first));
        
        _begin_pos = end+1;
        return (_begin_pos<=_end_pos);
    }

    void
    process_block(VcfRecord& vcfr,
                  const bool is_deblock,
                  const unsigned end) {

        if(! is_deblock) {
            if(end>vcfr.GetPos()) {
                vcfr.SetInfoVal("END",_intstr.get32(end));
            } else {
                vcfr.DeleteInfoKeyVal("END");
            }
            vcfr.WriteUnaltered(_opt.outfp);
        } else {
            vcfr.DeleteInfoKeyVal("END");
            vcfr.WriteUnaltered(_opt.outfp);
            while(end>vcfr.GetPos()) {
                const int next_pos(vcfr.GetPos()+1);
                vcfr.SetPos(next_pos);
                vcfr.SetRef(_scp.get_char(vcfr.GetChrom().c_str(),next_pos));
                vcfr.WriteUnaltered(_opt.outfp);
            }
        }
    }

    const BreakOptions& _opt;
    std::string _last_chrom;
    bool _is_skip_chrom; // true when pos is past all regions in current chrom
    region_util::interval_group_t::const_iterator _rhead,_rend;

    unsigned _begin_pos,_end_pos; // used to provide the region intercept iterator

    std::vector<int> _gti; // cache gt parse
    stringer<int> _intstr; // fast int->str util

    samtools_char_picker _scp;
};


static
void
process_vcf_input(const BreakOptions& opt,
                  std::istream& infp) {

    VcfHeaderHandler header(opt.outfp,gvcftools_version(),cmdline.c_str());
    BreakVcfRecordHandler rec(opt);

    istream_line_splitter vparse(infp);

    while(vparse.parse_line()) {
        if(header.process_line(vparse)) continue;
        rec.process_line(vparse);
    }
}






static
void
try_main(int argc,char* argv[]){

    //const time_t start_time(time(0));
    const char* progname(compat_basename(argv[0]));

    for(int i(0);i<argc;++i){
        if(i) cmdline += ' ';
        cmdline += argv[i];
    }

    std::istream& infp(std::cin);
    BreakOptions opt;
    std::string region_file;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
        ("region-file",po::value<std::string>(&region_file),"A bed file specifying regions where non-refernece blocks should be broken into individual positions (required)")
        ("ref", po::value<std::string >(&opt.ref_seq_file),"samtools reference sequence (required)");

    po::options_description help("help");
    help.add_options()
        ("help,h","print this message");

    po::options_description visible("options");
    visible.add(req).add(help);

    bool po_parse_fail(false);
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, visible), vm);
        po::notify(vm);    
    } catch(const boost::program_options::error& e) { // todo:: find out what is the more specific exception class thrown by program options
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }
    
    if ((argc<=1) || (vm.count("help")) || po_parse_fail) {
        log_os << "\n" << progname << " converts non-reference blocks to individual positions in specified regions\n\n"; 
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] < (g)VCF > unblocked_(g)VCF\n\n"; 
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }

    if(region_file.empty()) {
        log_os << "ERROR: no region file specified\n";
        exit(EXIT_FAILURE);
    }

    if(opt.ref_seq_file.empty()) {
        log_os << "ERROR: no region file specified\n";
        exit(EXIT_FAILURE);
    }

    region_util::get_regions(region_file,opt.regions);
    process_vcf_input(opt,infp);
}



static
void
dump_cl(int argc,
        char* argv[],
        std::ostream& os) {
 
    os << "cmdline:";
    for(int i(0);i<argc;++i){
        os << ' ' << argv[i];
    }
    os << std::endl;
}



int
main(int argc,char* argv[]){

    std::ios_base::sync_with_stdio(false);

    // last chance to catch exceptions...
    //
    try{
        try_main(argc,argv);

    } catch(const std::exception& e) {
        log_os << "FATAL:: EXCEPTION: " << e.what() << "\n"
               << "...caught in main()\n";
        dump_cl(argc,argv,log_os);
        exit(EXIT_FAILURE);

    } catch(...) {
        log_os << "FATAL:: UNKNOWN EXCEPTION\n"
               << "...caught in main()\n";
        dump_cl(argc,argv,log_os);
        exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
