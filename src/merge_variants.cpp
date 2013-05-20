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
    print_pos(const std::vector<boost::shared_ptr<site_crawler> >& sa,
              const reference_contig_segment& ref_seg);

private:

    std::ostream& _os;
    bool _is_header_output;
};



void
merge_reporter::
print_pos(const std::vector<boost::shared_ptr<site_crawler> >& sa,
          const reference_contig_segment& ref_seg)
{
    if(sa.empty()) return;

    if(! _is_header_output) {
        sa[0]->dump_header(_os);
        _os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\tFORMAT";

        const unsigned sample_size(sa.size());
        for(unsigned i(0);i<sample_size;++i) {
            _os << '\t' << sa[i]->sample_name();
        }
        _os << '\n';
        _is_header_output=true;
    }

    // start out getting ref, alt and gt for all samples:

    const char ref_base(ref_seg.get_base(sa[0]->pos-1));

    std::vector<std::string> genotypes;
    id_set<char> alleles;
    alleles.insert_key(ref_base);

    const unsigned n_samples(sa.size());
    for(unsigned st(0);st<n_samples;++st) {
        const site_crawler& sample(*(sa[st]));
        static const unsigned n_allele(2);
        std::ostringstream oss;
        for(unsigned ai(0); ai<n_allele; ai++) {
            const char allele(sample.allele[ai]);
            const unsigned code(alleles.insert_key(allele));
            oss << code;
            if(ai==0) oss << "/";
        }
        genotypes.push_back(oss.str());
    }

    std::ostringstream alt;
    const unsigned n_alleles(alleles.size());
    for(unsigned i(0);i<n_alleles;++i) {
        if(i) alt << ",";
        alt << alleles.get_key(i);
    }

    _os << sa[0]->chrom()
        << '\t' << sa[0]->pos
        << '\t' << '.' // ID
        << '\t' << ref_base
        << '\t' << alt.str()
        << '\t' << '.' // QUAL
        << '\t' << '.' // FILT
        << '\t' << '.' // INFO
        << '\t' << "GT"; // FORMAT

    const unsigned sample_size(genotypes.size());
    for(unsigned i(0);i<sample_size;++i) {
        _os << '\t' << genotypes[i];
    }

    _os << '\n';

#if 0
    for(unsigned i(0);i<sample_size;++i){
        *pos_fs_ptr << _sample_label[i] << "\t";
        *pos_fs_ptr << "\t";
        sa[i].dump_line(*pos_fs_ptr);
        *pos_fs_ptr << "\n";
    }
#endif
}


#if 0
struct site_stats : public site_stats_core<SAMPLE_SIZE> {

    site_stats()
    {
        for(unsigned i(0);i<PARENT_SIZE;++i) {
            for(unsigned j(0);j<CHILD_SIZE;++j) {
                snp_correct_type[i][j] = 0;
                incorrect_type[i][j] = 0;
            }
        }
    }

    unsigned snp_correct_type[PARENT_SIZE][CHILD_SIZE];
    unsigned incorrect_type[PARENT_SIZE][CHILD_SIZE];
};
#endif


static
void
merge_site(const std::vector<boost::shared_ptr<site_crawler> >& sa,
           const pos_t low_pos,
           const reference_contig_segment& ref_seg,
           merge_reporter& mr) {

    const unsigned n_samples(sa.size());

    bool is_all_mapped(true);
    bool is_any_mapped(false);

    bool is_all_called(true);
    bool is_any_called(false);
    
    bool is_any_nonref_called(false);

    const char ref_base(ref_seg.get_base(low_pos-1));

    if(ref_base=='N') return;

    for(unsigned st(0);st<n_samples;++st) {
        const site_crawler& site(*(sa[st]));
        if(site.is_pos_valid() && (site.pos==low_pos) && (site.n_total != 0)) {
            //ss.sample_mapped[st]++;
            is_any_mapped=true;
        } else {
            is_all_mapped=false;
        }

        if(site.is_pos_valid() && (site.pos==low_pos) && site.is_call){
            // position is called in sample st
            if(! ((ref_base==site.allele[0]) && (ref_base==site.allele[1]))){
                // position is called and non-ref at sample st
                is_any_nonref_called=true;
//                if(sa[st].allele[0] != sa[st].allele[1]){
//                    ss.sample_snp_het[st]++;
//                }
            }
            is_any_called=true;
        } else {
            is_all_called=false;
        }
    }

    // only interested in printing variants
    if(! is_any_nonref_called) return;

    
    mr.print_pos(sa,ref_seg);
}



#if 0
static
void
disallow_option(const boost::program_options::variables_map& vm,
                const char* label) {

    if(! vm.count(label)) return;

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
                             merge_reporter& mr) {

    // setup reference sequence:
    reference_contig_segment ref_seg;
    unsigned segment_known_size;
    get_samtools_std_ref_segment(ref_seq_file.c_str(),region,ref_seg,segment_known_size);
    if(opt.is_region()){
        ref_seg.set_offset(opt.region_begin-1);
    }

#if 0
    ss.ref_size += ref_seg.seq().size();
    ss.known_size += segment_known_size;
#endif

    // setup allele crawlers:
    std::vector<boost::shared_ptr<site_crawler> > sa;
    const unsigned n_samples(input_files.size());
    for(unsigned i(0);i<n_samples;++i){
        const bool is_store_header(i==0);
        sample_info tmp;
        tmp.file=input_files[i];
        sa.push_back(boost::shared_ptr<site_crawler>(new site_crawler(tmp, i, opt, region, ref_seg, is_store_header)));
    }

    while(true) {
        // get lowest position:
        bool is_low_pos_set(false);
        pos_t low_pos(0);
        for(unsigned st(0);st<n_samples;++st) {
            if(! sa[st]->is_pos_valid()) continue;
            if((! is_low_pos_set) || (sa[st]->pos < low_pos)){
                low_pos = sa[st]->pos;
                is_low_pos_set=true;
            }
        }
        if(! is_low_pos_set) break;
        
        merge_site(sa,low_pos,ref_seg,mr);

        for(unsigned st(0);st<n_samples;++st) {
            if(sa[st]->is_pos_valid() && (low_pos == sa[st]->pos)){
                sa[st]->update();
            }
        }
    }
}



static
void
try_main(int argc,char* argv[]){

//    const time_t start_time(time(0));
    const char* progname(compat_basename(argv[0]));
    
    for(int i(0);i<argc;++i){
        if(i) cmdline += ' ';
        cmdline += argv[i];
    }

    std::vector<std::string> input_files;
    std::string ref_seq_file;

    std::vector<std::string> exclude_list;

    shared_crawler_options opt;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
        ("ref", po::value<std::string>(&ref_seq_file),"samtools reference sequence (required)")
        ("region", po::value<std::string>(&opt.region), "samtools reference region (optional)")
        ("exclude", po::value<std::vector<std::string> >(&exclude_list), "name of chromsome to skip over (argument may be specified multiple times). Exclusions will be ignored if a region argument is provided")
        ("input", po::value<std::vector<std::string> >(&input_files)->multitoken(), "merge files (can be specified multiple times)");

    po::options_description help("help");
    help.add_options()
        ("help,h","print this message");

    po::options_description visible("options");
    visible.add(req).add(help);


    bool po_parse_fail(false);
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, visible,
                  po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);    
    } catch(const boost::program_options::error& e) { // todo:: find out what is the more specific exception class thrown by program options
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }
    
    bool is_show_help(false);
    if ((argc<=1) || (vm.count("help")) || po_parse_fail) {
        is_show_help=true;
    }

    if(ref_seq_file.empty()) is_show_help=true;
    if(input_files.empty()) is_show_help=true;

    if(is_show_help) {
        log_os << "\n" << progname << " merge the variants from multiple gVCF files\n\n";
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] > merged_variants\n\n";
        log_os << visible << "\n";
        exit(2);
    }

    if(opt.is_region()) {
        parse_tabix_region(input_files[0].c_str(),opt.region.c_str(),opt.region_begin,opt.region_end);
        opt.region_begin+=1;
    }

    merge_reporter mr(std::cout);
//    pos_reporters pr(conflict_pos_file,allhet_pos_file,hethethom_pos_file);
//    site_stats ss;

    if(opt.is_region()) {
        merge_variants(input_files,opt,ref_seq_file,opt.region.c_str(),mr);
    } else {
        fasta_chrom_list fcl(ref_seq_file.c_str());
        while(true) {
            const char* chrom = fcl.next();
            if(NULL == chrom) break;
            // don't even bother making this efficient:
            bool is_skip(false);
            for (unsigned i(0);i<exclude_list.size();++i) {
                if(strcmp(chrom,exclude_list[i].c_str())==0) {
                    is_skip=true;
                    break;
                }
            }
            if(is_skip) {
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