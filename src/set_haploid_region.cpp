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


struct SetHapOptions {

    SetHapOptions()
        : haploid_conflict_label("HAPLOID_CONFLICT")
        , orig_pl_tag("OPL")
        , outfp(std::cout)
    {}

    const std::string haploid_conflict_label;
    const std::string orig_pl_tag;
    std::ostream& outfp;

    region_util::region_t regions;
};



struct SetHapVcfHeaderHandler : public VcfHeaderHandler {

    SetHapVcfHeaderHandler(const SetHapOptions& opt,
                           const char* version = NULL,
                           const char* cmdline = NULL)
        : VcfHeaderHandler(opt.outfp,version,cmdline)
        , _opt(opt)
        , _is_add_filter_tag(true)
    {
        _haploid_filter_prefix="##FILTER=<ID="+opt.haploid_conflict_label;
    }

private:
    bool
    is_skip_header_line(const istream_line_splitter& vparse) {
        if(0 == strncmp(vparse.word[0],_haploid_filter_prefix.c_str(),_haploid_filter_prefix.size())) {
            _is_add_filter_tag=false;
        }
        return true;
    }

    void
    process_final_header_line() {
        if(_is_add_filter_tag) {
            _os << _haploid_filter_prefix
                << ",Description=\"Locus has heterozygous genotype in a haploid region.\">\n";
        }
        write_format(_opt.orig_pl_tag.c_str(),".","Integer","Original PL value before ploidy correction");
    }
    

    const SetHapOptions& _opt;
    bool _is_add_filter_tag;
    std::string _haploid_filter_prefix;
};



// process each vcf record for haploid setting:
//
struct SetHapVcfRecordHandler {

    SetHapVcfRecordHandler(const SetHapOptions& opt)
        : _opt(opt)
    {}

    void
    process_line(const istream_line_splitter& vparse) {
        const unsigned nw(vparse.n_word());

        if(nw <= VCFID::SAMPLE) {
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

        bool is_haploid;
        unsigned end;
        while(get_next_record_region_interval(is_haploid,end)){
            VcfRecord vcfr2(vcfr);
            if(end>vcfr2.GetPos()) {
                vcfr2.SetInfoVal("END",_intstr.get32(end));
            } else {
                vcfr2.DeleteInfoKeyVal("END");
            }
            if(is_haploid) make_record_haploid(vcfr2);
            vcfr2.WriteUnaltered(_opt.outfp);
            vcfr.SetPos(end+1);
        } 
        if(is_haploid) make_record_haploid(vcfr);
        vcfr.WriteUnaltered(_opt.outfp);
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
    get_next_record_region_interval(bool& is_haploid,
                                    unsigned& end) {

        assert(_begin_pos <= _end_pos);

        if(_begin_pos > _rhead->second) _rhead++;
        if(_rhead == _rend) { // no haploid regions left
            is_haploid = false;
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
        is_haploid = ((_begin_pos <= _rhead->second) && (end > _rhead->first));
        
        _begin_pos = end+1;
        return (_begin_pos<=_end_pos);
    }

    void
    make_record_haploid(VcfRecord& vcfr) {
        const char* gt(vcfr.GetSampleVal("GT"));
        if(NULL == gt)  return;
        parse_gt(gt,_gti);
            
        if(_gti.size() == 2) { // record is diploid
            if(_gti[0] == _gti[1]) {
                // change GT:
                static const char* unknown(".");
                const char* val(unknown);
                if(_gti[0]>=0) {
                    val=_intstr.get32(_gti[0]);
                }
                vcfr.SetSampleVal("GT",val);

                // move PL field to 'backup' OPL field:
                const char* pl(vcfr.GetSampleVal("PL"));
                if(NULL != pl) {
                    vcfr.SetSampleVal(_opt.orig_pl_tag.c_str(),pl);
                    vcfr.DeleteSampleKeyVal("PL");
                }
            } else {
                vcfr.AppendFilter(_opt.haploid_conflict_label.c_str());
            }
        }
    }


    const SetHapOptions& _opt;
    std::string _last_chrom;
    bool _is_skip_chrom; // true when pos is past all regions in current chrom
    region_util::interval_group_t::const_iterator _rhead,_rend;

    unsigned _begin_pos,_end_pos; // used to provide the region intercept iterator

    std::vector<int> _gti; // cache gt parse
    stringer<int> _intstr; // fast int->str util
};


static
void
process_vcf_input(const SetHapOptions& opt,
                  std::istream& infp) {

    SetHapVcfHeaderHandler header(opt,gvcftools_version(),cmdline.c_str());
    SetHapVcfRecordHandler rec(opt);

    istream_line_splitter vparse(infp);

    while(vparse.parse_line()) {
        if(header.process_line(vparse)) continue;

        if(vparse.n_word() > VCFID::SIZE) {
            std::ostringstream oss;
            oss << "Unexpected format in vcf record:\n";
            vparse.dump(oss);
            throw new blt_exception(oss.str().c_str());
        }

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
    SetHapOptions opt;
    std::string region_file;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
        ("region-file",po::value<std::string>(&region_file),"A bed file specifying the regions to be converted (required)");

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
        log_os << "\n" << progname << " converts regions of a gVCF or VCF from diploid to haploid\n\n"; 
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] < (g)VCF > haploid_region_(g)VCF\n\n"; 
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }

    if(region_file.empty()) {
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
