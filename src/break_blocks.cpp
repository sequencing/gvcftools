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

#include <unistd.h>

#include <iostream>
#include <string>


namespace {
std::ostream& log_os(std::cerr);
}

std::string cmdline;



// process each vcf record for haploid setting:
//
struct BreakVcfRecordHandler : public RegionVcfRecordHandler {

    BreakVcfRecordHandler(const RegionVcfOptions& opt)
        : RegionVcfRecordHandler(opt)
    {}

private:

    void
    process_block(const bool is_in_region,
                  const unsigned end,
                  VcfRecord& vcfr) const {

        if (! is_in_region) {
            if (end>vcfr.GetPos()) {
                vcfr.SetInfoVal("END",_intstr.get32(end));
            } else {
                vcfr.DeleteInfoKeyVal("END");
            }
            /// TODO: is it safe to pull the above if/else into this block?
            if (is_write_off_region_record(vcfr)) {
                vcfr.WriteUnaltered(_opt.outfp);
            }
        } else {
            vcfr.DeleteInfoKeyVal("END");
            vcfr.WriteUnaltered(_opt.outfp);
            while (end>vcfr.GetPos()) {
                const int next_pos(vcfr.GetPos()+1);
                vcfr.SetPos(next_pos);
                vcfr.SetRef(_scp.get_char(vcfr.GetChrom().c_str(),next_pos));
                vcfr.WriteUnaltered(_opt.outfp);
            }
        }
    }

    stringer<int> _intstr; // fast int->str util
};



static
void
process_vcf_input(const RegionVcfOptions& opt,
                  std::istream& infp) {

    VcfHeaderHandler header(opt.outfp,gvcftools_version(),cmdline.c_str());
    BreakVcfRecordHandler rec(opt);

    istream_line_splitter vparse(infp);

    while (vparse.parse_line()) {
        if (header.process_line(vparse)) continue;
        rec.process_line(vparse);
    }
}



static
void
try_main(int argc,char* argv[]) {

    //const time_t start_time(time(0));
    const char* progname(compat_basename(argv[0]));

    for (int i(0); i<argc; ++i) {
        if (i) cmdline += ' ';
        cmdline += argv[i];
    }

    std::istream& infp(std::cin);
    RegionVcfOptions opt;
    std::string region_file;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("region-file",po::value(&region_file),
     "A bed file specifying regions where call blocks should be broken into individual positions (required)")
    ("ref", po::value(&opt.refSeqFile),
     "samtools reference sequence (required)")
    ("exclude-off-target", po::value(&opt.isExcludeOffTarget)->zero_tokens(),
     "Don't output off-target vcf records. 'targeted' records include all those intersecting the input region plus any optionally included types specified below (default: output all records)")
    ("include-variants", po::value(&opt.isIncludeVariants)->zero_tokens(),
     "Add all variant calls to the targeted record set (only applies when exclude-off-target is used)");

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
    } catch (const boost::program_options::error& e) { // todo:: find out what is the more specific exception class thrown by program options
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }

    const bool isStdinTerminal(isatty(fileno(stdin)));

    if ((argc<=1) || (vm.count("help")) || po_parse_fail || isStdinTerminal) {
        log_os << "\n" << progname << " converts non-reference blocks to individual positions in specified regions\n\n";
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] < (g)VCF > unblocked_(g)VCF\n\n";
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }

    if (region_file.empty()) {
        log_os << "ERROR: no region file specified\n";
        exit(EXIT_FAILURE);
    }

    if (opt.refSeqFile.empty()) {
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
    for (int i(0); i<argc; ++i) {
        os << ' ' << argv[i];
    }
    os << std::endl;
}



int
main(int argc,char* argv[]) {

    std::ios_base::sync_with_stdio(false);

    // last chance to catch exceptions...
    //
    try {
        try_main(argc,argv);

    } catch (const std::exception& e) {
        log_os << "FATAL:: EXCEPTION: " << e.what() << "\n"
               << "...caught in main()\n";
        dump_cl(argc,argv,log_os);
        exit(EXIT_FAILURE);

    } catch (...) {
        log_os << "FATAL:: UNKNOWN EXCEPTION\n"
               << "...caught in main()\n";
        dump_cl(argc,argv,log_os);
        exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
