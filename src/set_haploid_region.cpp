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
#include "stringer.hh"
#include "VcfHeaderHandler.hh"
#include "VcfRecord.hh"


#include "boost/program_options.hpp"

//#include <ctime>

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>


namespace {
std::ostream& log_os(std::cerr);
}

std::string cmdline;


struct SetHapOptions {

    SetHapOptions()
        : haploid_conflict_label("HAPLOID_CONFLICT")
        , outfp(std::cout)
    {}


    const std::string haploid_conflict_label;
    std::ostream& outfp;

    typedef std::pair<unsigned,unsigned> interval_t;
    typedef std::vector<interval_t> interval_group_t;
    typedef std::map<std::string,interval_group_t> region_t;
    region_t regions;
};



static
void
parse_bedfile_regions(const std::string& region_file,
                      SetHapOptions::region_t& regions) {

    if (region_file.empty()) return;

    std::ifstream region_is(region_file.c_str());
    if (! region_is){
        log_os << "ERROR: Failed to open region file '" << region_file << "'\n";
        exit(EXIT_FAILURE);
    }

    unsigned line_no(0);
    bool is_parse_fail(false);

    std::string bed_chrom;
    unsigned bed_start(0),bed_end(0);

    while(! region_is.eof()){
        ++line_no;

        region_is >> bed_chrom;
        if(region_is.fail()) {
            if(! region_is.eof()) is_parse_fail=true;
            break;
        }

        if(bed_chrom != "track" && bed_chrom != "browser") {
        
            region_is >> bed_start >> bed_end;
            if(region_is.fail()) {
                is_parse_fail=true;
                break;
            }

            if(regions.find(bed_chrom) == regions.end()) {
                regions[bed_chrom] = SetHapOptions::interval_group_t();
            }
            regions[bed_chrom].push_back(std::make_pair(bed_start,bed_end));
        }

        region_is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if(is_parse_fail) {
        log_os << "ERROR: unexpected format in region bed file line no: " << line_no << "\n";
        exit(EXIT_FAILURE);
    }
}



/// sort and merge overlapping/adjacent regions:
static
void
merge_regions(SetHapOptions::region_t& regions) {

    SetHapOptions::interval_group_t tmp_regions;

    SetHapOptions::region_t::iterator chromi(regions.begin()), chromi_end(regions.end());
    for(;chromi!=chromi_end;++chromi) {
        SetHapOptions::interval_group_t& cregions(chromi->second);
        std::sort(cregions.begin(),cregions.end());

        const unsigned nr(cregions.size());
        for(unsigned i(0);i<nr;++i) {
            SetHapOptions::interval_t& ci(cregions[i]);
            SetHapOptions::interval_t& hi(i==0 ? ci : tmp_regions.back());
            if((i==0) || (ci.first > hi.second)) {
                tmp_regions.push_back(ci);
            } else if(ci.second > hi.second) {
                hi = std::make_pair(hi.first,ci.second);
            }
        }
        cregions.swap(tmp_regions);
        tmp_regions.clear();
    }
}



void
dump_regions(const SetHapOptions::region_t& regions,
             std::ostream& os) {

    SetHapOptions::region_t::const_iterator chromi(regions.begin()), chromi_end(regions.end());
    for(;chromi!=chromi_end;++chromi) {
        const std::string& chrom(chromi->first);
        const SetHapOptions::interval_group_t& cregions(chromi->second);

        const unsigned nr(cregions.size());
        for(unsigned i(0);i<nr;++i) {
            os << chrom << "\t" << cregions[i].first << "\t" << cregions[i].second << "\n";
        }
    }
}



static
void
get_regions(const std::string& region_file,
            SetHapOptions::region_t& regions) {

    parse_bedfile_regions(region_file,regions);
    merge_regions(regions);
}



struct SetHapVcfHeaderHandler : public VcfHeaderHandler {

    SetHapVcfHeaderHandler(const SetHapOptions& opt,
                           const char* version = NULL,
                           const char* cmdline = NULL)
        : VcfHeaderHandler(opt.outfp,version,cmdline)
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
    }
    

    bool _is_add_filter_tag;
    std::string _haploid_filter_prefix;
};



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

        if(! is_in_region(vparse)) {
            vparse.write_line(_opt.outfp);
            return;
        }

        // write modified record:
        VcfRecord vcfr(vparse);    
        const char* gt(vcfr.GetSampleVal("GT"));
        if(NULL != gt) {
            parse_gt(gt,_gti);
            
            if(_gti.size() == 2) {
                if(_gti[0] == _gti[1]) {
                    vcfr.SetSampleVal("GT",_stringer.itos_32(_gti[0]));
                } else {
                    vcfr.AppendFilter(_opt.haploid_conflict_label.c_str());
                }
            }
        }
        vcfr.WriteUnaltered(_opt.outfp);
   }

private:
    bool
    is_in_region(const istream_line_splitter& vparse) {
        // determine if chromosome is new:
        if(_last_chrom.empty() || (0 != strcmp(_last_chrom.c_str(),vparse.word[0]))) {
            _last_chrom=vparse.word[0];
            const SetHapOptions::region_t::const_iterator i(_opt.regions.find(_last_chrom));
            _is_skip_chrom=((i==_opt.regions.end()) || (i->second.empty()));

            // setup region iterators:
            if(! _is_skip_chrom) {
                _rhead=i->second.begin();
                _rend=i->second.end();
            }
        }

        if(! _is_skip_chrom) {
            const char* posstr(vparse.word[1]);
            const unsigned pos(parse_unsigned(posstr));
            while(_rhead != _rend) {
                if(pos>_rhead->second) {
                    _rhead++;
                } else {
                    return (pos>_rhead->first);
                }
            }
            _is_skip_chrom=true;
        }
        return false;
    }

    const SetHapOptions& _opt;
    std::string _last_chrom;
    bool _is_skip_chrom; // true when pos is past all regions in current chrom
    SetHapOptions::interval_group_t::const_iterator _rhead,_rend;

    std::vector<int> _gti; // cache gt parse
    stringer _stringer; // fast int->str util
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
        log_os << "\n" << progname << " convert regions of a gVCF or VCF from diploid to haploid\n\n"; 
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] < all_sites > trio_report\n\n"; 
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }

    if(region_file.empty()) {
        log_os << "ERROR: no region file specified\n";
        exit(EXIT_FAILURE);
    }

    get_regions(region_file,opt.regions);
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
