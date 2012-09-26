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
#include "RegionVcfRecordHandler.hh"
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


struct SetHapOptions : public RegionVcfOptions {

    SetHapOptions()
        : haploid_conflict_label("HAPLOID_CONFLICT")
        , orig_pl_tag("OPL")
    {}

    const std::string haploid_conflict_label;
    const std::string orig_pl_tag;
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
struct SetHapVcfRecordHandler : public RegionVcfRecordHandler {

    SetHapVcfRecordHandler(const SetHapOptions& opt)
        : RegionVcfRecordHandler(opt)
        , _shopt(opt)
    {}

private:

    void
    process_block(const bool is_in_region,
                  const unsigned end,
                  VcfRecord& vcfr) const {

        if(end>vcfr.GetPos()) {
            vcfr.SetInfoVal("END",_intstr.get32(end));
        } else {
            vcfr.DeleteInfoKeyVal("END");
        }
        if(is_in_region) make_record_haploid(vcfr);
        vcfr.WriteUnaltered(_opt.outfp);
    }

    void
    make_record_haploid(VcfRecord& vcfr) const {
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
                    vcfr.SetSampleVal(_shopt.orig_pl_tag.c_str(),pl);
                    vcfr.DeleteSampleKeyVal("PL");
                }
            } else {
                vcfr.AppendFilter(_shopt.haploid_conflict_label.c_str());
            }
        }
    }

    const SetHapOptions& _shopt;
    mutable std::vector<int> _gti; // cache gt parse
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
        ("region-file",po::value<std::string>(&region_file),"A bed file specifying the regions to be converted (required)")
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

    if(opt.ref_seq_file.empty()) {
        log_os << "ERROR: no reference file specified\n";
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
