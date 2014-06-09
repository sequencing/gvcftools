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

#include "gvcftools.hh"
#include "id_map.hh"
#include "related_sample_util.hh"
#include "ref_util.hh"
#include "string_util.hh"
#include "tabix_util.hh"
#include "trio_option_util.hh"

#include "boost/foreach.hpp"
#include "boost/program_options.hpp"
#include "boost/shared_ptr.hpp"

#include <ctime>

#include <iostream>
#include <string>
#include <vector>


std::ostream& log_os(std::cerr);
std::ostream& report_os(std::cout);

std::string cmdline;



struct merge_reporter {

    merge_reporter(std::ostream& os)
        : _os(os),
          _is_header_output(false)
    {}

    void
    print_locus(
        const std::vector<boost::shared_ptr<site_crawler> >& sa,
        const vcf_pos low_pos,
        const reference_contig_segment& ref_seg,
        const bool is_indel=false);


    static const char filter_delim;
    static const char format_delim;

private:

    std::ostream& _os;
    bool _is_header_output;

    // not persistent, just used to reduce allocation:
    std::vector<std::string> words;
};


const char merge_reporter::filter_delim(';');
const char merge_reporter::format_delim(':');



template <typename T>
void
refAltWriter(
    const id_set<T>& alleles,
    std::ostream& os)
{
    const unsigned n_alleles(alleles.size());
    assert(0 != n_alleles);

    os << '\t' << alleles.get_key(0);  // REF

    // ALT:
    os << '\t';
    if (n_alleles>1) {
        for (unsigned i(1); i<n_alleles; ++i) {
            if (i>1) os << ",";
            os << alleles.get_key(i);
        }
    } else {
        os << '.';
    }
}



void
merge_reporter::
print_locus(
    const std::vector<boost::shared_ptr<site_crawler> >& sa,
    const vcf_pos low_pos,
    const reference_contig_segment& ref_seg,
    const bool is_indel)
{
    if (sa.empty()) return;

    if (! _is_header_output) {
        sa[0]->dump_header(_os);
        _os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\tFORMAT";

        const unsigned sample_size(sa.size());
        for (unsigned i(0); i<sample_size; ++i) {
            _os << '\t' << sa[i]->sample_name();
        }
        _os << '\n';
        _is_header_output=true;
    }

    const unsigned n_samples(sa.size());

    //
    // merge all filters
    //
    typedef id_set<std::string> fkey_t;
    fkey_t merged_filters;
    for (unsigned st(0); st<n_samples; ++st) {
        const site_crawler& sample(*(sa[st]));
        if (!(low_pos == sample.vpos())) continue;
        const char* filter(sample.word(VCFID::FILT));
        if ((0 != strcmp(filter,".")) && (0 != strcmp(filter,"PASS"))) {
            split_string(filter,filter_delim,words);
            BOOST_FOREACH(const std::string& w, words)
            {
                merged_filters.insert_key(w);
            }
        }
    }

    //
    // merge all format fields
    //
    fkey_t merged_keys;
    std::vector<fkey_t> sample_keys(n_samples);
    for (unsigned st(0); st<n_samples; ++st) {
        const site_crawler& sample(*(sa[st]));
        if (!(low_pos == sample.vpos())) continue;
        const char* format(sample.word(VCFID::FORMAT));
        if (0 != strcmp(format,".")) {
            split_string(format,format_delim,words);
            BOOST_FOREACH(const std::string& w, words)
            {
                merged_keys.insert_key(w);
                sample_keys[st].insert_key(w);
            }
        }
    }

    _os << sa[0]->chrom()      // CHROM
        << '\t' << low_pos.pos // POS
        << '\t' << '.';        // ID


    std::vector<std::string> genotypes;

    if (! is_indel)
    {
        // get ref, alt and gt for all samples:
        //
        const char ref_base(ref_seg.get_base(sa[0]->pos()-1));

        id_set<char> alleles;
        alleles.insert_key(ref_base);

        for (unsigned st(0); st<n_samples; ++st) {
            const site_crawler& sample(*(sa[st]));
            bool is_nonstandard(false);
            std::ostringstream oss;
            if (low_pos == sample.vpos()) {
                const unsigned n_allele(sample.get_allele_size());

                for (unsigned ai(0); ai<n_allele; ai++) {
                    const char allele(sample.get_allele(ai));
                    if (allele != 'N') {
                        const unsigned code(alleles.insert_key(allele));
                        if (ai) oss << '/';
                        oss << code;
                    } else {
                        is_nonstandard=true;
                    }
                }
            }
            else
            {
                is_nonstandard=true;
            }
            if (is_nonstandard) {
                genotypes.push_back(".");
            } else {
                genotypes.push_back(oss.str());
            }
        }

        refAltWriter(alleles,_os);
    }
    else
    {   //indel case:

        // get ref, alt and gt for all samples:
        //
        std::string ref_allele;
        for (unsigned st(0); st<n_samples; ++st) {
            const site_crawler& sample(*(sa[st]));
            if (!(low_pos == sample.vpos())) continue;
            if (sample.get_indel_ref().size() > ref_allele.size())
            {
                ref_allele=sample.get_indel_ref();
            }
        }

        std::vector<std::string> alt_adjust(n_samples);
        for (unsigned st(0); st<n_samples; ++st) {
            const site_crawler& sample(*(sa[st]));
            if (!(low_pos == sample.vpos())) continue;
            if (sample.get_indel_ref().size() < ref_allele.size())
            {
                alt_adjust[st] = ref_allele.substr(sample.get_indel_ref().size());
            }
        }

        id_set<std::string> alleles;
        alleles.insert_key(ref_allele);

        for (unsigned st(0); st<n_samples; ++st) {
            const site_crawler& sample(*(sa[st]));
            if (low_pos == sample.vpos()) {
                const unsigned n_allele(sample.get_indel_allele_size());

                std::ostringstream oss;
                for (unsigned ai(0); ai<n_allele; ai++) {
                    const std::string allele(sample.get_indel_allele(ai));
                    if (allele != "X") {
                        const unsigned code(alleles.insert_key(allele+alt_adjust[st]));
                        if (ai) oss << '/';
                        oss << code;
                    } else {
                        if (ai) oss << '/';
                        oss << 0;
                    }
                }
                genotypes.push_back(oss.str());
            }
            else
            {
                genotypes.push_back(".");
            }
        }

        refAltWriter(alleles,_os);
    }

    assert(genotypes.size() == n_samples);


    _os << '\t' << '.'; // QUAL

    // FILT:
    _os << '\t';
    const unsigned n_filters(merged_filters.size());
    if (n_filters) {
        for (unsigned filter_index(0); filter_index<n_filters; ++filter_index) {
            if (filter_index) _os << filter_delim;
            _os << merged_filters.get_key(filter_index);
        }
    } else {
        _os << "PASS";
    }

    _os << '\t' << '.'; // INFO

    // FORMAT:
    _os << '\t';
    const unsigned n_keys(merged_keys.size());
    if (n_keys) {
        for (unsigned i(0); i<n_keys; ++i) {
            if (i) _os << format_delim;
            _os << merged_keys.get_key(i);
        }
    } else {
        _os << '.';
    }

    for (unsigned st(0); st<n_samples; ++st) {
        const site_crawler& sample(*(sa[st]));
        const fkey_t& sample_key(sample_keys[st]);

        if (low_pos == sample.vpos()) {
            const char* vcf_sample(sample.word(VCFID::SAMPLE));

            if (0 != strcmp(vcf_sample,".")) {
                split_string(vcf_sample,format_delim,words);
            } else {
                words.clear();
            }
        }
        else
        {
            words.clear();
        }

        // print out values in merged order:
        _os << '\t';
        if (n_keys && (low_pos == sample.vpos())) {
            for (unsigned key_index(0); key_index<n_keys; ++key_index) {
                if (key_index) _os << format_delim;
                const std::string& key(merged_keys.get_key(key_index));
                if (sample_key.test_key(key)) {
                    if (key == "GT") {
                        _os << genotypes[st];
                    } else {
                        _os << words[sample_key.get_id(key)];
                    }
                } else {
                    _os << '.';
                }
            }
        } else {
            _os << '.';
        }
    }

    _os << '\n';

#if 0
    for (unsigned i(0); i<sample_size; ++i) {
        *pos_fs_ptr << _sample_label[i] << "\t";
        *pos_fs_ptr << "\t";
        sa[i].dump_line(*pos_fs_ptr);
        *pos_fs_ptr << "\n";
    }
#endif
}



static
void
merge_site(const std::vector<boost::shared_ptr<site_crawler> >& sa,
           const vcf_pos low_pos,
           const reference_contig_segment& ref_seg,
           merge_reporter& mr) {

    if (low_pos.is_indel) return;

    const unsigned n_samples(sa.size());

    bool is_any_nonref_called(false);

    const char ref_base(ref_seg.get_base(low_pos.pos-1));

    if (ref_base=='N') return;

    for (unsigned st(0); st<n_samples; ++st) {
        const site_crawler& site(*(sa[st]));

        if (site.is_pos_valid() && (site.vpos()==low_pos) && site.is_site_call()) {
            // position is called in sample st
            const unsigned n_allele(site.get_allele_size());
            for (unsigned allele_index(0); allele_index<n_allele; allele_index++) {
                if (ref_base != site.get_allele(allele_index)) {
                    is_any_nonref_called=true;
                }
            }
        }
    }

    // only interested in printing variants
    if (! is_any_nonref_called) return;

    mr.print_locus(sa, low_pos, ref_seg);
}



#if 0
static
void
disallow_option(const boost::program_options::variables_map& vm,
                const char* label) {

    if (! vm.count(label)) return;

    log_os << "ERROR:: option '" << label << "' is not allowed in selected snp-mode\n";
    exit(EXIT_FAILURE);
}
#endif




static
void
merge_variants(const std::vector<std::string>& input_files,
               const shared_crawler_options& opt,
               const std::string& ref_seq_file,
               const char* region,
               merge_reporter& mr)
{
    // setup reference sequence:
    reference_contig_segment ref_seg;
    unsigned segment_known_size;
    get_samtools_std_ref_segment(ref_seq_file.c_str(),region,ref_seg,segment_known_size);
    if (opt.is_region()) {
        ref_seg.set_offset(opt.region_begin-1);
    }

#if 0
    ss.ref_size += ref_seg.seq().size();
    ss.known_size += segment_known_size;
#endif

    // setup locus crawlers:
    std::vector<boost::shared_ptr<site_crawler> > sa;
    const unsigned n_samples(input_files.size());
    for (unsigned i(0); i<n_samples; ++i) {
        static const bool is_return_indels(true);
        const bool is_store_header(i==0);
        sample_info tmp;
        tmp.file=input_files[i];
        sa.push_back(boost::shared_ptr<site_crawler>(new site_crawler(tmp, i, opt, region, ref_seg, is_store_header, is_return_indels)));
    }

    while (true) {
        // get lowest position:
        bool is_low_pos_set(false);
        vcf_pos low_pos;
        for (unsigned st(0); st<n_samples; ++st) {
            if (! sa[st]->is_pos_valid()) continue;
            if ((! is_low_pos_set) || (sa[st]->vpos() < low_pos)) {
                low_pos = sa[st]->vpos();
                is_low_pos_set=true;
            }
        }
        if (! is_low_pos_set) break;

        if (! low_pos.is_indel) {
            merge_site(sa, low_pos, ref_seg, mr);
        }
        else
        {
            mr.print_locus(sa, low_pos, ref_seg,true);
        }

        for (unsigned st(0); st<n_samples; ++st) {
            if (sa[st]->is_pos_valid() && (low_pos == sa[st]->vpos())) {
                sa[st]->update();
            }
        }
    }
}



static
void
try_main(int argc,char* argv[]) {

//    const time_t start_time(time(0));
    const char* progname(compat_basename(argv[0]));

    for (int i(0); i<argc; ++i) {
        if (i) cmdline += ' ';
        cmdline += argv[i];
    }

    std::vector<std::string> input_files;
    std::string ref_seq_file;

    std::vector<std::string> exclude_list;

    shared_crawler_options opt;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("ref", po::value(&ref_seq_file),"samtools reference sequence (required)")
    ("region", po::value(&opt.region), "samtools reference region (optional)")
    ("exclude", po::value<std::vector<std::string> >(&exclude_list), "name of chromosome to skip over (argument may be specified multiple times). Exclusions will be ignored if a region argument is provided")
    ("input", po::value<std::vector<std::string> >(&input_files)->multitoken(), "merge files (can be specified multiple times)")
    ("murdock", po::value(&opt.is_murdock_mode)->zero_tokens(),
     "If true, don't stop because of any out-of-order position conflicts. Any out of order positions are ignored. In case of an overlap the first observation is used and subsequent repeats are ignored.")
    ;

    po::options_description help("help");
    help.add_options()
    ("help,h","print this message");

    po::options_description visible("options");
    visible.add(req).add(help);


    bool po_parse_fail(false);
    po::variables_map vm;
    try {
        po::parsed_options parsed(po::parse_command_line(argc, argv, visible,
                                  po::command_line_style::unix_style ^ po::command_line_style::allow_short));
        po::store(parsed,vm);

        po::notify(vm);

        // any remaining options are an error:
        if (! po::collect_unrecognized(parsed.options,po::include_positional).empty())
        {
            log_os << "\nERROR: Unexpected positional options.\n\n";
            po_parse_fail=true;
        }
    } catch (const boost::program_options::error& e) { // todo:: find out what is the more specific exception class thrown by program options
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n\n";
        po_parse_fail=true;
    }

    bool is_show_help(false);
    if ((argc<=1) || (vm.count("help")) || po_parse_fail) {
        is_show_help=true;
    }

    if (ref_seq_file.empty()) is_show_help=true;
    if (input_files.empty()) is_show_help=true;

    if (is_show_help) {
        log_os << "\n" << progname << " merges the variants from multiple gVCF files\n\n";
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] > merged_variants\n\n";
        log_os << visible << "\n";
        exit(2);
    }

    if (opt.is_region()) {
        parse_tabix_region(input_files[0].c_str(),opt.region.c_str(),opt.region_begin,opt.region_end);
        opt.region_begin+=1;
    }

    merge_reporter mr(std::cout);
//    pos_reporters pr(conflict_pos_file,allhet_pos_file,hethethom_pos_file);
//    site_stats ss;

    if (opt.is_region()) {
        merge_variants(input_files,opt,ref_seq_file,opt.region.c_str(),mr);
    } else {
        fasta_chrom_list fcl(ref_seq_file.c_str());
        while (true) {
            const char* chrom = fcl.next();
            if (NULL == chrom) break;
            // don't even bother making this efficient:
            bool is_skip(false);
            for (unsigned i(0); i<exclude_list.size(); ++i) {
                if (strcmp(chrom,exclude_list[i].c_str())==0) {
                    is_skip=true;
                    break;
                }
            }
            if (is_skip) {
                log_os << "skipping chromosome: '" << chrom << "'\n";
            } else {
                log_os << "processing chromosome: '" << chrom << "'\n";
                merge_variants(input_files,opt,ref_seq_file,chrom,mr);
            }
        }

    }
//    report(ss,start_time,is_variable_metadata);
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
