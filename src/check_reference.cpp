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
#include "stringer.hh"
#include "VcfHeaderHandler.hh"
#include "vcf_util.hh"

#include "boost/program_options.hpp"

#include <unistd.h>

#include <iostream>
#include <string>


namespace {
std::ostream& log_os(std::cerr);
}

std::string cmdline;


struct RefCheckOptions
{
    RefCheckOptions() : outfp(std::cout) {}

    std::ostream& outfp;
    std::string refSeqFile;
};



// process each vcf record for haploid setting:
//
struct RefCheckVcfRecordHandler
{
    RefCheckVcfRecordHandler(const RefCheckOptions& opt)
        : _opt(opt)
    {}

    void
    process_line(const istream_line_splitter& vparse)
    {
        const unsigned nw(vparse.n_word());

        if (nw < (VCFID::INFO+1)) {
            log_os << "ERROR: unexpected number of fields in vcf record:\n";
            vparse.dump(log_os);
            exit(EXIT_FAILURE);
        }

        const char* chrom = vparse.word[VCFID::CHROM];

        const char* pos_ptr(vparse.word[VCFID::POS]);
        const unsigned pos = parse_unsigned(pos_ptr);
        assert(pos>0);

        const char* ref = vparse.word[VCFID::REF];
        if (*ref == '\0')
        {
            std::ostringstream oss;
            oss << "Empty reference field in input vcf record:\n";
            vparse.dump(oss);
            throw blt_exception(oss.str().c_str());
        }

        if (_last_chrom.empty() || (0 != strcmp(chrom,_last_chrom.c_str())))
        {
            unsigned known_size(0);
            get_samtools_std_ref_segment(_opt.refSeqFile.c_str(),chrom,_ref, known_size);
            _last_chrom = chrom;
        }

        const unsigned refSize(strlen(ref));
        const std::string& refstr(_ref.seq());
        if (0 != refstr.compare(pos-1,refSize,ref))
        {
            std::cerr << "ERROR: vcf REF value '" << ref << "' conflicts with fasta ref value '" << refstr.substr(pos-1,refSize) <<"'. At vcf line:\n";
            vparse.dump(std::cerr);
            exit(EXIT_FAILURE);
        }
    }


private:
    const RefCheckOptions& _opt;
    std::string _last_chrom;
    reference_contig_segment _ref;
};



static
void
process_vcf_input(
    const RefCheckOptions& opt,
    std::istream& infp)
{
    VcfHeaderHandler header(opt.outfp,gvcftools_version(),cmdline.c_str(),true);
    RefCheckVcfRecordHandler rec(opt);

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
    RefCheckOptions opt;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("ref", po::value(&opt.refSeqFile),
     "samtools reference sequence (required)")
    ;

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
        log_os << "\n" << progname << " check VCF reference fields.\n\n";
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] < (g)VCF\n\n";
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }

    if (opt.refSeqFile.empty()) {
        log_os << "ERROR: no reference file specified\n";
        exit(EXIT_FAILURE);
    }

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
