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
#include "related_sample_util.hh"
#include "ref_util.hh"
#include "tabix_util.hh"
#include "trio_option_util.hh"

#include "boost/program_options.hpp"

#include <ctime>

#include <iostream>
#include <string>
#include <vector>

std::ostream& log_os(std::cerr);
std::ostream& report_os(std::cout);

std::string cmdline;


enum sample_t {
    TWIN1,
    TWIN2,
    SAMPLE_SIZE
};

const char* sample_label[] = { "twin1","twin2" };


enum state_t {
    SAMEHOM,
    DIFFHOM,
    SAMEHET,
    DIFFHET,
    HOMHET,
    STATE_SIZE
};

const char* state_label[] = { "samehom" , "diffhom" , "samehet" , "diffhet" , "homhet" };


static
state_t
get_st_cat(const bool ist1hom,
           const bool ist2hom,
           const char t1_1,
           const char t2_1,
           const char t1_2,
           const char t2_2){

    if(ist1hom!=ist2hom) return HOMHET;

    if(ist1hom) {
        return (t1_1==t2_1 ? SAMEHOM : DIFFHOM);
    } else {
        return ((t1_1==t2_1) && (t1_2==t2_2) ? SAMEHET : DIFFHET);
    }
}



struct site_stats : public site_stats_core<SAMPLE_SIZE> {

    site_stats() {
        for(unsigned i(0);i<STATE_SIZE;++i) {
            snp_correct_type[i]=0;
            incorrect_type[i]=0;
        }
    }

    unsigned snp_correct_type[STATE_SIZE];
    unsigned incorrect_type[STATE_SIZE];
};



static
void
processSite(const site_crawler* sa,
            const pos_t low_pos,
            const reference_contig_segment& ref_seg,
            pos_reporter& pr,
            site_stats& ss) {

    bool is_all_mapped(true);
    bool is_any_mapped(false);

    bool is_all_called(true);
    bool is_any_called(false);
    
    const char ref_base(ref_seg.get_base(low_pos-1));

    if(ref_base=='N') return;

    for(unsigned st(0);st<SAMPLE_SIZE;++st) {
        if(sa[st].is_pos_valid() && (sa[st].pos()==low_pos) && (sa[st].n_total() != 0)) {
            ss.sample_mapped[st]++;
            is_any_mapped=true;
        } else {
            is_all_mapped=false;
        }

        if(sa[st].is_pos_valid() && (sa[st].pos()==low_pos) && sa[st].is_call()){
            ss.sample_called[st]++;
            if(! ((ref_base==sa[st].get_allele(0)) && (ref_base==sa[st].get_allele(1)))){
                ss.sample_snp[st]++;
                if(sa[st].get_allele(0) != sa[st].get_allele(1)){
                    ss.sample_snp_het[st]++;
                }
            }
            is_any_called=true;
        } else {
            is_all_called=false;
        }
    }

    if(! is_all_mapped) {
        if(is_any_mapped) ss.some_mapped++;
    } else {
        ss.all_mapped++;
    }

    if(! is_all_called) {
        if(is_any_called) ss.some_called++;
        return;
    } else {
        ss.all_called++;
    }

    const char t1_1(sa[TWIN1].get_allele(0));
    const char t1_2(sa[TWIN1].get_allele(1));
    const char t2_1(sa[TWIN2].get_allele(0));
    const char t2_2(sa[TWIN2].get_allele(1));

    const bool is_correct(((t1_1==t2_1) && (t1_2==t2_2)) ||
                          ((t1_2==t2_1) && (t1_1==t2_2)));

    const bool ist1hom(t1_1==t1_2);
    const bool ist2hom(t2_1==t2_2);


    if(is_correct) {
        const bool is_ref_call(ist1hom && ist2hom && (t1_1==ref_base));
        if(! is_ref_call){

            const state_t st(get_st_cat(ist1hom,ist2hom,t1_1,t2_1,t1_2,t2_2));

            ss.snp_correct++;
            ss.snp_correct_type[st]++;

            if     (! ist1hom) ss.sample_snp_correct_het[TWIN1]++;
            else if(t1_1!=ref_base) ss.sample_snp_correct_hom[TWIN1]++;
            if     (! ist2hom) ss.sample_snp_correct_het[TWIN2]++;
            else if(t2_1!=ref_base) ss.sample_snp_correct_hom[TWIN2]++;
        }
    } else {
        pr.print_pos(sa);

        const state_t st(get_st_cat(ist1hom,ist2hom,t1_1,t2_1,t1_2,t2_2));

        ss.incorrect++;
        ss.incorrect_type[st]++;
    }
}



static
void
report(const site_stats& ss,
       const time_t& start_time,
       const bool is_variable_metadata){

    std::ostream& os(report_os);

    os << "CMDLINE " << cmdline << "\n";
    if(is_variable_metadata) {
        os << "START_TIME " << asctime(localtime(&start_time));
        os << "VERSION " << gvcftools_version() << "\n";
    }
    os << "\n";
    os << "sites: " << ss.ref_size << "\n";
    os << "known_sites: " << ss.known_size << "\n";
    os << "\n";
    for(unsigned st(0);st<SAMPLE_SIZE;++st) {
        os << "sites_mapped_" << sample_label[st] << ": " << ss.sample_mapped[st] << "\n";
    }
    assert(ss.known_size >= (ss.some_mapped+ss.all_mapped));
    const unsigned none_mapped(ss.known_size-(ss.some_mapped+ss.all_mapped));
    os << "sites_mapped_in_no_samples: " << none_mapped << "\n";
    os << "sites_mapped_in_some_samples: " << ss.some_mapped << "\n";
    os << "sites_mapped_in_all_samples: " << ss.all_mapped << "\n";
    os << "\n";
    for(unsigned st(0);st<SAMPLE_SIZE;++st) {
        os << "sites_called_" << sample_label[st] << ": " << ss.sample_called[st] << "\n";
    }
    assert(ss.known_size >= (ss.some_called+ss.all_called));
    const unsigned none_called(ss.known_size-(ss.some_called+ss.all_called));
    os << "sites_called_in_no_samples: " << none_called << "\n";
    os << "sites_called_in_some_samples: " << ss.some_called << "\n";
    os << "sites_called_in_all_samples: " << ss.all_called << "\n";
    os << "sites_called_in_all_samples_conflict: " << ss.incorrect << "\n";
    os << "fraction_of_known_sites_called_in_all_samples: " << ratio(ss.all_called,ss.known_size) << "\n";
    os << "fraction_of_sites_called_in_all_samples_in_conflict: " << ratio(ss.incorrect,ss.all_called) << "\n";
    os << "\n";
    const unsigned snps(ss.snp_correct+ss.incorrect);
    for(unsigned st(0);st<SAMPLE_SIZE;++st) {
        const unsigned het(ss.sample_snp_het[st]);
        const unsigned hom(ss.sample_snp[st]-ss.sample_snp_het[st]);
        const double sample_het_hom(ratio(het,hom));  
        const double sample_phet(ratio(het,(hom+het)));  
        os << "sites_with_snps_called_total_het_hom_het/hom_P(het)_" << sample_label[st] << ": " << ss.sample_snp[st] << " " << het << " " << hom << " " << sample_het_hom << " " << sample_phet << "\n";
    }
    os << "sites_called_in_all_samples_with_snps_called_any_sample: " << snps << "\n";
    os << "fraction_of_snp_sites_in_conflict: " << ratio(ss.incorrect,snps) << "\n";
    os << "\n";

    for(unsigned st(0);st<SAMPLE_SIZE;++st) {
        const unsigned het(ss.sample_snp_correct_het[st]);
        const unsigned hom(ss.sample_snp_correct_hom[st]);
        const double sample_het_hom(ratio(het,hom));  
        const double sample_phet(ratio(het,(hom+het)));  
        os << "snp_non_conflict_total_het_hom_het/hom_P(het)_" << sample_label[st] << ": " << (het+hom) << " " << het << " " << hom << " " << sample_het_hom << " " << sample_phet << "\n";
    }
    os << "\n";

    for(unsigned i(0);i<STATE_SIZE;++i) {
        os << "snp_conflict_type_" << state_label[i] << ": "
               << ss.incorrect_type[i] << "\n";
    }
    os << "\n";
    for(unsigned i(0);i<STATE_SIZE;++i) {
        os << "snp_non_conflict_type_" << state_label[i] << ": "
           << ss.snp_correct_type[i] << "\n";
    }
    os << "\n";
}



static
void
accumulate_region_statistics(const sample_info* const si,
                             const shared_crawler_options& opt,
                             const std::string& ref_seq_file,
                             const char* region,
                             pos_reporter& pr,
                             site_stats& ss) {

    // setup reference sequence:
    reference_contig_segment ref_seg;
    unsigned segment_known_size;
    get_samtools_std_ref_segment(ref_seq_file.c_str(),region,ref_seg,segment_known_size);
    if(opt.is_region()){
        ref_seg.set_offset(opt.region_begin-1);
    }

    ss.ref_size += ref_seg.seq().size();
    ss.known_size += segment_known_size;

    // setup allele crawlers:
    site_crawler sa[SAMPLE_SIZE] = { site_crawler(si[TWIN1],TWIN1,opt,region,ref_seg),
                                     site_crawler(si[TWIN2],TWIN2,opt,region,ref_seg) };

    while(true) {
        // get lowest position:
        bool is_low_pos_set(false);
        pos_t low_pos(0);
        for(unsigned st(0);st<SAMPLE_SIZE;++st) {
            if(! sa[st].is_pos_valid()) continue;
            if((! is_low_pos_set) || (sa[st].pos() < low_pos)){
                low_pos = sa[st].pos();
                is_low_pos_set=true;
            }
        }
        if(! is_low_pos_set) break;

        processSite(sa,low_pos,ref_seg,pr,ss);

        for(unsigned st(0);st<SAMPLE_SIZE;++st) {
            if(sa[st].is_pos_valid() && (low_pos == sa[st].pos())){
                sa[st].update();
            }
        }
    }
}



static
void
try_main(int argc,char* argv[]){

    const time_t start_time(time(0));
    const char* progname(compat_basename(argv[0]));
    
    for(int i(0);i<argc;++i){
        if(i) cmdline += ' ';
        cmdline += argv[i];
    }

    std::string conflict_pos_file;
    std::string ref_seq_file;
    shared_crawler_options opt;

    std::vector<std::string> exclude_list;

    sample_info si[SAMPLE_SIZE];
    snp_param& sp(opt.sp());

    std::vector<info_filter> max_info_filters;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
        ("ref", po::value<std::string >(&ref_seq_file),"samtools reference sequence (required)")
        ("region", po::value<std::string>(&opt.region), "samtools reference region (optional)")
        ("exclude", po::value<std::vector<std::string> >(&exclude_list), "name of chromsome to skip over (argument may be specified multiple times). Exclusions will be ignored if a region argument is provided")
        ("twin1", po::value<std::string>(&si[TWIN1].file), 
         "twin/replicate 1 gvcf file")
        ("twin2", po::value<std::string>(&si[TWIN2].file), 
         "twin/replicate 2 gvcf file")
        ("conflict-file", po::value<std::string>(&conflict_pos_file), "Write all conflict positions to the specified file")
        ("no-variable-metadata",
         "Remove timestamp and any other metadata from output during validation testing")
        ("murdock",
         "In murdock mode twins doesn't stop because of any out-of-order position conflicts. Any out of order positions are ignored. In case of an overlap the first observation is used and subsequent repeats are ignored.");

    po::options_description filter("filtration");
    filter.add_options()
        ("min-gqx", po::value<double>(&sp.min_gqx), "If GQX value for a record is below this value, then don't use the locus. Note that if the filter field already contains a GQX filter, this will not 'rescue' filtered variants when min-gqx is set very low -- this filter can only lower callability on a file. Any records missing the GQX field will not be filtered out. (default: 0)")
        ("min-pos-rank-sum", po::value<double>(&sp.min_pos_rank_sum), "Filter site if the INFO field contains the key BaseQRankSum and the value is less than the minimum. (default: no-filter)")
        ("min-qd", po::value<double>(&sp.min_qd), "Filter site if the INFO field contains the key QD and the value is less than the minimum. (default: no-filter)")
        ("min-info-field",po::value<std::vector<info_filter> >(&sp.infof)->multitoken(),
         "Filter records which contain an INFO key equal to argument1, and a corresponding value less than argument2 ")
        ("max-info-field",po::value<std::vector<info_filter> >(&max_info_filters)->multitoken(),
         "Filter records which contain an INFO key equal to argument1, and a corresponding value greater than argument2 ");

    po::options_description help("help");
    help.add_options()
        ("help,h","print this message");

    po::options_description visible("options");
    visible.add(req).add(filter).add(help);

    bool po_parse_fail(false);
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, visible,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);    
    } catch(const boost::program_options::error& e) { // todo:: find out what is the more specific exception class thrown by program options
        log_os << "ERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }
    
    if ((argc<=1) || (vm.count("help")) || ref_seq_file.empty() || po_parse_fail) {
        log_os << "\n" << progname << " finds conflicts in the variant calls made from twins or technical replicates.\n\n"; 
        log_os << "version: " << gvcftools_version() << "\n\n";
        log_os << "usage: " << progname << " [options] > twins_report\n\n"; 
        log_os << visible << "\n";
	log_os << "Note that calls inside of deletions will not be used\n";
        exit(EXIT_FAILURE);
    }

    // clean up filters:
    {
        for(unsigned i(0);i<max_info_filters.size();++i) {
            max_info_filters[i].is_min=false;
            sp.infof.push_back(max_info_filters[i]);
        }
    }

    const bool is_variable_metadata(! vm.count("no-variable-metadata"));
    opt.is_murdock_mode = vm.count("murdock");

    sp.is_min_qd=(vm.count("min-qd"));
    sp.is_min_pos_rank_sum=(vm.count("min-pos-rank-sum"));

#if 0
    for(unsigned st(0);st<SAMPLE_SIZE;++st) {
        log_os << sample_label[st] << "_files:\n";
        const unsigned st_size(si[st].allele_files.size());
        for(unsigned i(0);i<st_size;++i) {
            log_os << si[st].allele_files[i] << "\n";
        }
    }
#endif

    for(unsigned st(0);st<SAMPLE_SIZE;++st) {
        if(si[st].file.empty()) {
            log_os << "ERROR: no gvcf file specified for sample: '" << sample_label[st] << "'\n";
            exit(EXIT_FAILURE);
        }
    }

    if(opt.is_region()) {
        parse_tabix_region(si[TWIN1].file.c_str(),opt.region.c_str(),opt.region_begin,opt.region_end);
        opt.region_begin+=1;
    }

    std::vector<std::string> slabel(sample_label,sample_label+SAMPLE_SIZE);
    pos_reporter pr(conflict_pos_file,slabel);
    site_stats ss;

    if(opt.is_region()) {
        accumulate_region_statistics(si,opt,ref_seq_file,opt.region.c_str(),pr,ss);
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
                accumulate_region_statistics(si,opt,ref_seq_file,chrom,pr,ss);
            }
        }

    }
    report(ss,start_time,is_variable_metadata);
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
