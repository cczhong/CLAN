# include "mapping_driver.h"

using namespace std;

void MappingDriver::ReadMapBatch(
    char **ref_header, char **ref_seq, int num_refs,
    char **header, char **seq, int num_reads,          // the header and the sequence of the read set
    const int begin, const int size, const int upper,  // the begin, the number of sequences, and the upper bound of sequences to map in this batch
    BWT &bwt,                           // the indexed BWT
    std::vector< std::vector<FragType> > &result,      // the output
    const CLANSearchParam &param_set
)
{
    assert(size > 0);
    //cerr << "begin mapping" << endl;
////////////////////////////////////////
    //cerr << "num refs: " << num_refs << endl;
    //BWTIDX *ref_len = new BWTIDX [num_refs];
    //for(int k = 0; k < num_refs; ++ k) {
    //  string s = ref_seq[k];
    //  ref_len[k] = (BWTIDX) s.length();
    //}
    //delete [] ref_len;
    //return;
////////////////////////////////////////


    for(int i = 0; i < size && begin + i < upper; ++ i)
    {
        //cerr << "Working on " << header[begin + i] << endl;
        //cerr << "seq  " << seq[begin + i] << endl;
        BWTSearch searcher;
        BWTFragMerger merger;

        vector<AlignType> all_hits, all_hits_rv;
        vector<MatchType> merged_frag, extended_frag, merged_frag_rv, extended_frag_rv;
        string s = seq[begin + i];

        // search forward direction
        searcher.SearchAllSubRegions(bwt, param_set.min_frag_len, param_set.max_hits, param_set.flank_len, s.c_str(), all_hits);
        //continue;
        searcher.ExtendAllSubRegions(bwt, s.c_str(), all_hits, extended_frag);
        //continue;
        merger.MergeFragments(bwt, s.length(), extended_frag, merged_frag);
        //continue;

        //cerr << "frag_size: " << extended_frag.size() << endl;

        for(auto it = merged_frag.begin(); it != merged_frag.end(); ++ it)
        {
            it->strand = true;
        }

        //continue;

        //cerr << s << endl;
        //cerr << "reverse search begins" << endl;
        // if the data is not strand-specific, map to both strands
        if(!param_set.strand_specific)
        {
            string s_rv = string(s.rbegin(), s.rend());
            // search reverse direction
            for(int j = 0; j < s_rv.length(); ++ j)
            {
                switch(s_rv[j])
                {
                case 'A':
                    s_rv[j] = 'T';
                    break;
                case 'C':
                    s_rv[j] = 'G';
                    break;
                case 'G':
                    s_rv[j] = 'C';
                    break;
                case 'N':
                    s_rv[j] = 'N';
                    break;
                case 'T':
                    s_rv[j] = 'A';
                    break;
                default:
                    break;
                }
            }
            //cerr << s_rv << endl;
            searcher.SearchAllSubRegions(bwt, param_set.min_frag_len, param_set.max_hits, param_set.flank_len, s_rv.c_str(), all_hits_rv);
            searcher.ExtendAllSubRegions(bwt, s_rv.c_str(), all_hits_rv, extended_frag_rv);
            merger.MergeFragments(bwt, s_rv.length(), extended_frag_rv, merged_frag_rv);
            // add the newly detected fragments into the candidate pool
            for(auto it = merged_frag_rv.begin(); it != merged_frag_rv.end(); ++ it)
            {
                it->strand = false;
                int b = s.length() - it->q_end - 1;
                int e = s.length() - it->q_begin - 1;
                it->q_begin = b;
                it->q_end = e;
                merged_frag.push_back(*it);
            }
        }
        //cerr << "before best matching" << endl;
        //continue;
        int covered = merger.FindBestMatching(
                          bwt, param_set.num_frag, merged_frag, result[i],
                          param_set.use_insert_penalty, param_set.frag_penalty, param_set.max_overlap
                      );

        //cerr << "i: " << i << endl;

////////////////////////////////////////
//    for(int k = 0; k < result[i].size(); ++ k) {
//      for(int l = 0; l < result[i][k].targets.size(); ++ l) {
//        if(result[i][k].targets[l].sid >= num_refs)  {
//          cerr << "ref ID error in search!!!" << result[i][k].targets[l].sid << "  " << num_refs << endl;
//        }
//        if(result[i][k].targets[l].end >= ref_len[result[i][k].targets[l].sid])  {
//          cerr << "Search error!!!" << header[begin + i] << endl;
//        }
//      }
//    }
////////////////////////////////////////
    }


    return;
}


