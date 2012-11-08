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

#include "BlockerOptions.hh"
#include "BlockerVcfHeaderHandler.hh"
#include "blt_exception.hh"
#include "gvcftools.hh"
#include "istream_line_splitter.hh"
#include "parse_util.hh"
#include "VcfRecordBlocker.hh"

#include "boost/program_options.hpp"

//#include <ctime>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


namespace {
std::ostream& log_os(std::cerr);
}

std::string cmdline;



static
void
process_vcf_input(const BlockerOptions& opt,
                  std::istream& infp) {

    VcfRecordBlocker blocker(opt);
    BlockerVcfHeaderHandler header(opt,gvcftools_version(),cmdline.c_str());

    istream_line_splitter vparse(infp);

    while(vparse.parse_line()) {
        if(header.process_line(vparse)) continue;

        if(vparse.n_word() > VCFID::SIZE) {
            std::ostringstream oss;
            oss << "Unexpected format in vcf record:\n";
            vparse.dump(oss);
            throw new blt_exception(oss.str().c_str());
        }

        try {
            GatkVcfRecord record(vparse);
            blocker.Append(record);
        } catch (const std::exception& e) {
            log_os << "ERROR: Exception thrown while processing vcf record: '" << e.what() << "'\n"
                   << "\tVCF_INPUT_STATE:\n";
            vparse.dump(log_os);
            log_os << "\n";
            throw;
        }
    }
}



// parse the chrom depth file
static
void
parse_chrom_depth(const std::string& chrom_depth_file,
                  std::map<std::string, double>& ChromDepth) {

    if (chrom_depth_file.empty()) return;

    std::ifstream depth_is(chrom_depth_file.c_str());
    if (! depth_is){
        log_os << "ERROR: Failed to open chrom depth file '" << chrom_depth_file << "'\n";
        exit(EXIT_FAILURE);
    }

    static const unsigned buff_size(1024);
    char buff[buff_size];

    unsigned line_no(0);

    while(true){
        depth_is.getline(buff,buff_size);
        if(! depth_is) {
            if     (depth_is.eof()) break;
            else {
                log_os << "ERROR: unexpected failure while attempting to read chrom depth file line " << (line_no+1) << "\n";
                exit(EXIT_FAILURE);
            }
        } else {
            ++line_no;
        }

        char* word2(strchr(buff,'\t'));
        if(NULL == word2) {
            log_os << "ERROR: unexpected format in read chrom depth file line " << (line_no) << "\n";
            exit(EXIT_FAILURE);
        }
        *(word2++) = '\0';
        try {
            const char* s(word2);
            ChromDepth[buff] = parse_double(s);
        } catch(const blt_exception& e) {
            log_os << "ERROR: unexpected format in read chrom depth file line " << (line_no) << "\n";
            throw;
        }
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
    BlockerOptions opt;
    std::string chrom_depth_file;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
        ("min-blockable-nonref",po::value<print_double>(&opt.min_nonref_blockable)->default_value(opt.min_nonref_blockable),"If AD present, only compress non-variant site if 1-AD[0]/DP < value")
        ("skip-header","Write gVCF output without header");

    po::options_description filters("filters");
    filters.add_options()
        ("chrom-depth-file",po::value<std::string>(&chrom_depth_file),"Read mean depth for each chromosome from file, and use these values for maximum site depth filteration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)")
        ("max-depth-factor",po::value<print_double>(&opt.max_chrom_depth_filter_factor)->default_value(opt.max_chrom_depth_filter_factor),"If a chrom depth file is supplied then loci with depth exceeding the mean chrom depth times this value are filtered")
        ("min-gqx",po::value<std::string>(&opt.min_gqx)->default_value(opt.min_gqx),"Minimum locus GQX");

    for(unsigned i(0);i<opt.filters.size();++i) {
        FilterInfo& fi(opt.filters[i]);
        filters.add_options()
            (fi.argname.c_str(),po::value<print_double>(&fi.thresh)->default_value(fi.thresh),fi.GetArgDescription().c_str());
    }

    filters.add_options()
        ("no-default-filters","Clear all default filters. Any individual filter threshold changes above will still be in effect");

    po::options_description blocks("blocks");
    blocks.add_options()
        ("block-range-factor",po::value<print_double>(&opt.nvopt.BlockFracTol)->default_value(opt.nvopt.BlockFracTol),
         "Non-variant blocks are restricted to range [x,y], y <= max(x+3,x*(1+block-range-factor))")
        ("block-label",po::value<std::string>(&opt.nvopt.BlockavgLabel)->default_value(opt.nvopt.BlockavgLabel),
         "VCF INFO key used to annotate compressed non-variant blocks");

    po::options_description help("help");
    help.add_options()
        ("help,h","print this message");

    po::options_description visible("options");
    visible.add(filters).add(req).add(blocks).add(help);

    bool po_parse_fail(false);
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, visible), vm);
        po::notify(vm);    
    } catch(const boost::program_options::error& e) { // todo:: find out what is the more specific exception class thrown by program options
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }

    //    if ((argc<=1) || (vm.count("help")) || po_parse_fail) {
    if ((vm.count("help")) || po_parse_fail) {
        log_os << "\n" << progname << " creates block-compressed gVCF from modified GATK all sites output\n\n"; 
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] < all_sites > gVCF\n\n"; 
        log_os << visible << "\n";
        exit(2);
    }

    opt.is_skip_header=vm.count("skip-header");

    if(opt.nvopt.BlockFracTol.numval() < 0) {
        log_os << "\nblock-range-factor must be >= 0\n\n";
        exit(2);
    }

    if(vm.count("no-default-filters")) {
        if(vm["min-gqx"].defaulted()) opt.min_gqx.clear();

        std::vector<FilterInfo> new_filters;
        for(unsigned i(0);i<opt.filters.size();++i) {
            const FilterInfo& fi(opt.filters[i]);
            if(! vm[fi.argname.c_str()].defaulted()) new_filters.push_back(fi);
        }
        opt.filters = new_filters;
    }

    if(! chrom_depth_file.empty()) {
        parse_chrom_depth(chrom_depth_file,opt.ChromDepth);
    }

    opt.finalize_filters();

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
