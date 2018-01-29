#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <ctime>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <omp.h>

#include "bwt.h"
#include "bwt_search.h"
#include "loader.h"
#include "bio_alphabet.h"
#include "bwt_frag_merger.h"
#include "mapping_driver.h"
#include "output_printer.h"
#include "parameters.h"
#include "misc.h"

#ifndef NUM_FRAGS
#define NUM_FRAGS 2
#endif

#ifndef MAP_CHUNK_SIZE
#define MAP_CHUNK_SIZE 1000
#endif

using namespace std;


//int num_threads = 1;
//int max_hits = 20;
//int min_frag_len = 10;
//int flank_len = 4;
//int num_frag = 2;
//int frag_penalty = 10;
//int max_overlap = 5;
//bool strand_specific = false;
//bool use_insert_penalty = true;
//char *f_reference = new char [1000];
//char *f_index = new char [1000];
//char *f_input = new char [1000];
//char *f_output = new char [1000];

static CLANSearchParam param_set;
int main(int argc, char **argv)
{

    // setting default parameters
    param_set.num_threads = 1;
    param_set.max_hits = 20;
    param_set.min_frag_len = 10;
    param_set.flank_len = 4;
    param_set.num_frag = 2;
    param_set.frag_penalty = 10;
    param_set.max_overlap = 5;
    param_set.strand_specific = false;
    param_set.use_insert_penalty = true;
    param_set.f_reference[0] = '\0';
    param_set.f_index[0] = '\0';
    param_set.f_input[0] = '\0';
    param_set.f_output[0] = '\0';

    int copt;
    extern char *optarg;
    extern int optind;
    while ((copt=getopt(argc,argv,"r:o:f:d:t:m:n:k:l:p:v:seh")) != EOF)
    {
        switch(copt)
        {
        case 'r':
            sscanf(optarg, "%s", param_set.f_input);
            continue;
        case 'o':
            sscanf(optarg, "%s", param_set.f_output);
            continue;
        case 'f':
            sscanf(optarg, "%s", param_set.f_reference);
            continue;
        case 'd':
            sscanf(optarg, "%s", param_set.f_index);
            continue;
        case 't':
            sscanf(optarg, "%d", &param_set.num_threads);
            continue;
        case 'm':
            sscanf(optarg, "%d", &param_set.max_hits);
            continue;
        case 'k':
            sscanf(optarg, "%d", &param_set.flank_len);
            continue;
        case 'l':
            sscanf(optarg, "%d", &param_set.min_frag_len);
            continue;
        case 'p':
            sscanf(optarg, "%d", &param_set.frag_penalty);
            continue;
        case 'v':
            sscanf(optarg, "%d", &param_set.max_overlap);
            continue;
        case 's':
            param_set.strand_specific = true;
            continue;
        case 'e':
            param_set.use_insert_penalty = false;
            continue;
        case 'h':
        default:
            cout << "==========================================================" << endl;
            cout << "\tCLAN: the CrossLinked reads ANalaysis tool" << endl;
            cout << "==========================================================" << endl;
            cout << endl;
            cout << "usage: clan_search -r [READ_FILE] -o [OUTPUT_FILE] -f [REFERENCE_FILE] -d [INDEX_PREFIX]" << endl;
            cout << endl;
            cout << "\tr: the file containing the reads, in FASTA format (mandatory)" << endl;
            cout << "\to: the file for writing the temporary mapping results (mandatory)" << endl;
            cout << "\tf: the reference sequence (e.g. the human genome, mandatory)" << endl;
            cout << "\td: the prefix of the indexes (files built by \"clan_index\", mandatory)" << endl;
            cout << "\ts: enable strand-specific mapping (ignore mapping to reverse strand, default FALSE)" << endl;
            cout << "\te: disable insert penalty (penalize insert sequence between two arms, default TRUE)" << endl;
            cout << "\tt: number of threads to use (optional, default 1)" << endl;
            cout << "\tm: number of maximum hits for each maximal fragment (optional, default 20)" << endl;
            cout << "\tk: length of flanking bases when recording non-maximal fragments (optional, default 4)" << endl;
            cout << "\tl: minimum length for each fragment (optional, default 10)" << endl;
            cout << "\tp: penalty for introducing one extra strand (optional, default 10)" << endl;
            cout << "\tv: maximum overlap allowed between mapped regions (optional, default 5)" << endl;
            cout << "\th: print this help message" << endl << endl;
            exit(0);
        }
        optind--;
    }

    if(strlen(param_set.f_reference) <= 0 || strlen(param_set.f_index) <= 0 ||
            strlen(param_set.f_input) <= 0 || strlen(param_set.f_output) <= 0)
    {
        cerr << "Mandatory argument missing; please type \"clan_search -h\" to view the help information." << endl;
        cerr << "Abort." << endl;
        exit(1);
    }
    ofstream out_fh;
    out_fh.open(param_set.f_output, ios::out);
    if(!out_fh.is_open())
    {
        cerr << "clan_search: unable to write to " << param_set.f_output << "; Abort." << endl;
        exit(1);
    }

    //out_fh << "#read_id\tsolution_id\tread_mapped_begin\tread_mapped_end\tread_length\tmapped_locations" << endl;

    double start_time = misc::MyTime();
    double check_time;

    BioAlphabet dna_alphabet(DNA);
    Loader loader;

    int num_refs = loader.CountFastaNumSeqs(param_set.f_reference);
    char **ref_header = new char* [num_refs];
    char **ref_seq = new char* [num_refs];
    num_refs = loader.LoadFasta(dna_alphabet, param_set.f_reference, ref_header, ref_seq);
    cerr << "CLAN: Finish loading reference" << endl;

    string concat_seq;
    Concatenator concat_obj(ref_seq, num_refs, concat_seq);

    //for(int i = 0; i < concat_seq.length() - 13; ++ i) {
    //  if(concat_seq.substr(i, 13) == "GGGTTCGAGCCCC") {
    //    cerr << "found occurrence" << endl;
    //  }
    //}

    BWT bwt;
    bwt.ConstructFromIndex(
        dna_alphabet, concat_seq.c_str(), param_set.f_index
    );
    //cout << "Done constructing object" << endl;
    bwt.ConstructLocationInfo(num_refs, ref_header, ref_seq);
    cerr << "CLAN: Finish constructing BWT" << endl;


    // load sequece data

    int num_reads;
    if(loader.IsFASTA(param_set.f_input))
    {
        num_reads = loader.CountFastaNumSeqs(param_set.f_input);
    }
    else if(loader.IsFASTQ(param_set.f_input))
    {
        num_reads = loader.CountFastqNumSeqs(param_set.f_input);
    }
    else
    {
        cerr << "CLAN: unrecognized reads format, only FASTA or FASTQ format is suported. Abort." << endl;
        exit(1);
    }

    cerr << "Num reads: " << num_reads << endl;

    char **header = new char* [num_reads];
    char **seq = new char* [num_reads];
    if(loader.IsFASTA(param_set.f_input))
    {
        num_reads = loader.LoadFasta(dna_alphabet, param_set.f_input, header, seq);
    }
    else
    {
        num_reads = loader.LoadFastq(dna_alphabet, param_set.f_input, header, seq);
    }

    check_time = misc::MyTime();
    //PrintElapsed(start_time, check_time, "Loading index");
    start_time = misc::MyTime();
    cerr << "CLAN: Finish loading reads" << endl;


    BWTIDX num_hits = 0;
    // perform the search

    vector<vector<FragType> > selected_frag;
    vector<int> seq_order;

    MappingDriver map_driver;
    OutputPrinter out_printer;

    int num_chunks = 0;
    for(int i = 0; i < num_reads; i += MAP_CHUNK_SIZE)
    {
        ++ num_chunks;
    }
    out_fh.write((char*) &num_chunks, sizeof(int));


    #pragma omp parallel num_threads(param_set.num_threads)
    {
        #pragma omp for
        for(int i = 0; i < num_reads; i += MAP_CHUNK_SIZE)
        {

            vector<vector<FragType> > *mapping_result = new vector<vector<FragType> >;
            mapping_result->resize(MAP_CHUNK_SIZE);
            // clearing up the vector
            for(int j = 0; j < MAP_CHUNK_SIZE; ++ j)
            {
                (*mapping_result)[j].clear();
            }
            //cerr << "Entering mapping" << endl;
            map_driver.ReadMapBatch(
                ref_header, ref_seq, num_refs,
                header, seq, num_reads,
                i, MAP_CHUNK_SIZE, num_reads,
                bwt, *mapping_result,
                param_set
            );
            # pragma omp critical
            {
                //cerr << "Preparing for output" << endl;
                out_printer.OutputEncoded(
                    ref_seq, num_refs,
                    i, MAP_CHUNK_SIZE, num_reads,
                    *mapping_result, out_fh
                );
                //cerr << "Output done" << endl;
                delete mapping_result;
            }
        }
    }
    out_fh.close();
    check_time = misc::MyTime();
    //PrintElapsed(start_time, check_time, "Writing all results");

    //cerr << "Preparing for input" << endl;
    //vector<vector<FragType> > check_load;
    //out_printer.ReadEncoded(
    //  ref_seq, num_refs,
    //  num_reads, string(param_set.f_output), check_load
    //);
    //cerr << "Input done" << endl;

    for(int i = 0; i < num_refs; ++ i)
    {
        delete [] ref_header[i];
        delete [] ref_seq[i];
    }
    delete [] ref_header;
    delete [] ref_seq;
    for(int i = 0; i < num_reads; ++ i)
    {
        delete [] header[i];
        delete [] seq[i];
    }
    delete [] header;
    delete [] seq;

    return 0;
}
