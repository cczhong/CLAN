#include "bwt_frag_merger.h"

using namespace std;



bool SortTargetTypeByBegin(const TargetType &a, const TargetType &b)
{
    if(a.begin < b.begin || (a.begin == b.begin && a.end < b.end))
    {
        return true;
    }
    return false;
}

bool SortMatchTypeByQBegin(const MatchType &a, const MatchType &b)
{
    if(a.q_begin < b.q_begin || (a.q_begin == b.q_begin && a.q_end > b.q_end) ||
            (a.q_begin == b.q_begin && a.q_end == b.q_end && a.str_pos < b.str_pos)
      )
    {
        return true;
    }
    return false;
}

bool SortMatchType(const MatchType &a, const MatchType &b)
{
    if(a.str_pos < b.str_pos || (a.str_pos == b.str_pos && a.q_begin < b.q_begin) ||
            (a.str_pos == b.str_pos && a.q_begin == b.q_begin && a.q_end > b.q_end)
      )
    {
        return true;
    }
    return false;
}

bool SortAlignTypeByBWTPos(const AlignType &a, const AlignType &b)
{
    if(a.bwt_begin < b.bwt_begin || (a.bwt_begin == b.bwt_begin && a.bwt_end < b.bwt_end))
    {
        return true;
    }
    return false;
}

/*
void BWTFragMerger::AlignToFragType(
    BWT &bwt, std::vector<AlignType> &all_pos, std::vector<FragType> &all_frag, const int num_hits
) {
  if(all_pos.size() <= 0) return;
  for(auto it = all_pos.begin(); it != all_pos.end(); ++ it) {
    FragType f; f.q_begin = it->q_begin; f.q_end = it->q_end;
    long long int diff = (long long int) (it->bwt_end) - (long long int) (it->bwt_begin);
    //cout << "!!!!!  " << diff << " " << num_hits << " " << (long long int) it->bwt_end << " " << (long long int) it->bwt_begin << endl;
    if(diff > num_hits) {
      //cout << "skipped..." << endl;
      continue;
    }
    //cout << "region: " << f.q_begin << " " << f.q_end << endl;
    //cout << "index: " << it->bwt_begin << " " << it->bwt_end << endl;
    for(BWTIDX i = it->bwt_begin; i < it->bwt_end; ++ i) {
      //cout << "died at " << i << endl;
      TargetType t;
      bwt.GetRefLocation(i, t.sid, t.begin);
      t.end = t.begin + f.q_end - f.q_begin;
      cout << f.q_begin << "  " << f.q_end << " " << t.sid << "  " << t.begin << " " << t.end << endl;
      f.targets.push_back(t);
    }
    all_frag.push_back(f);
    //cout << "hello world" << endl;
  }
  //cout << "num processed fragments: " << all_frag.size() << endl;
  return;
}
*/

void BWTFragMerger::MergeFragments(
    BWT &bwt, const int query_len,
    std::vector<MatchType> &frags, std::vector<MatchType> &merged_frags,
    const int gap
)
{
    // sort the matchings based on the positions in the reference sequences
    sort(frags.begin(), frags.end(), SortMatchType);

    //for(auto it = frags.begin(); it != frags.end(); ++ it) {
    //cerr << (long long int) it->str_pos << "  " << it->q_begin << " " << it->q_end << endl;
    //  cerr << it->cigar << endl;
    //}

    // parse the matchings and define the fragments
    bool *q_pos = new bool [DIFF_ARRAY_SIZE];
    BWTIDX last = -100;
    BWTIDX chunk_begin = 2;
    vector<MatchType> nr_match;
    for(int i = 0; i < frags.size(); ++ i)
    {
        if(frags[i].str_pos != last && frags[i].str_pos - last != 1)
        {
            memset(q_pos, false, DIFF_ARRAY_SIZE);
            chunk_begin = frags[i].str_pos;
        }
        last = frags[i].str_pos;
        int diff_idx = (frags[i].str_pos - frags[i].q_begin) - (chunk_begin - query_len);
        if(diff_idx < 0 || diff_idx >= DIFF_ARRAY_SIZE)
        {
            cerr << "Warning: BWTFragMerger::MergeFragments: diff array index out of range!!!" << endl;
            continue;
        }
        if(!q_pos[diff_idx])
        {
            q_pos[diff_idx] = true;
            //cout << "!!!  " << (long long int) frags[i].str_pos << "  " << frags[i].q_begin << " " << frags[i].q_end << endl;
            nr_match.push_back(frags[i]);
        }
        else
        {
            delete [] frags[i].cigar;
        }
    }
    delete [] q_pos;

    //for(auto it = nr_match.begin(); it != nr_match.end(); ++ it) {
    //  string ref_id; BWTIDX ref_idx;
    //  bwt.StrPosToRefPos(it->str_pos, ref_id, ref_idx);
    //  cout << ":::  " << (long long int) it->str_pos << "  " << ref_id << "  " << (long long int) ref_idx << " " << it->q_begin << " " << it->q_end << endl;
    //  cerr << "???: " << it->q_begin << " " << it->q_end << endl;
    //}

    // attemp to progressively merge hits
    // NOTE: we assume that the hits in nr_match are sorted

    for(int i = 0; i < nr_match.size(); ++ i)
    {
        if(nr_match[i].str_len <= 0) continue;
        int last_q_end = nr_match[i].q_end;
        BWTIDX right_bound = nr_match[i].str_pos + nr_match[i].str_len - 1;
        int acc_score = nr_match[i].score;
        int acc_mismatch = nr_match[i].num_mismatch;
        char *merged_cigar = new char [50];
        strcpy(merged_cigar, nr_match[i].cigar);
        //cerr << "i: " << nr_match[i].q_begin << " " << nr_match[i].q_end << endl;
        // record the last record
        if(i == nr_match.size() - 1)
        {
            merged_frags.push_back(nr_match[i]);
            break;
        }
        for(int j = i + 1; j < nr_match.size(); ++ j)
        {
            //cerr << "j: " << nr_match[j].q_begin << " " << nr_match[j].q_end << endl;
            // if this the two blocks can be merged
            int g_str = nr_match[j].str_pos - right_bound - 1;
            int g_qry = nr_match[j].q_begin - last_q_end - 1;
            int gap_size = abs(g_str - g_qry);
            char operation = g_str < g_qry ? 'I' : 'D';
            int mismatch_size = g_str < g_qry ? g_str : g_qry;
            int gap_score = mismatch_size * NT_MISMATCH + (gap_size > 0 ? NT_GAP_OPEN + gap_size * NT_GAP_EXT : 0);
            //cerr << "criteria:  " << nr_match[j].str_len << " " << mismatch_size << " " << gap_size << endl;
            if(nr_match[j].str_len > 0 && mismatch_size <= gap && gap_size <= gap
                    && nr_match[j].score + gap_score > 0
                    && g_str >= 0 && g_qry >= 0
              )
            {
                //cerr << "merged" << endl;
                // record merged block

                if(mismatch_size > 0) sprintf(merged_cigar, "%s%d%c", merged_cigar, mismatch_size, 'M');
                if(gap_size > 0) sprintf(merged_cigar, "%s%d%c", merged_cigar, gap_size, operation);
                sprintf(merged_cigar, "%s%s", merged_cigar, nr_match[j].cigar);

                last_q_end = nr_match[j].q_end;
                right_bound = nr_match[j].str_pos + nr_match[j].str_len - 1;
                //cerr << "score check: " << acc_score << " " << nr_match[j].score << " " << gap_score << endl;
                acc_score += nr_match[j].score + gap_score;
                acc_mismatch += nr_match[j].num_mismatch;
                nr_match[j].str_len = 0;
                delete [] nr_match[j].cigar;
            }
            else if(g_str > gap || g_qry > gap)
            {
                break;
            }
        }
        // record the current matching
        MatchType t;
        t.str_pos = nr_match[i].str_pos;
        t.str_len = right_bound - t.str_pos + 1;
        t.q_begin = nr_match[i].q_begin;
        t.q_end = last_q_end;
        t.score = acc_score;
        t.num_mismatch = acc_mismatch;
        t.cigar = new char [50];
        strcpy(t.cigar, merged_cigar);

        //t.cigar = merged_cigar;
        //cerr << merged_cigar << " " << t.q_begin << "  " << t.q_end << " " << t.str_pos << " " << t.str_len << endl;
        merged_frags.push_back(t);
        delete [] merged_cigar;
        delete [] nr_match[i].cigar;
    }

    return;
}


int BWTFragMerger::FindBestMatching(
    BWT &bwt, const int num_fragments,
    std::vector<MatchType> &all_frag,
    std::vector<FragType> &selected,
    const bool use_insert_penalty,
    const int frag_cost, const int max_overlap
)
{
    // filter out matchings with too many hits
    int i, j, k;
    sort(all_frag.begin(), all_frag.end(), SortMatchTypeByQBegin);

    //for(auto it = all_frag.begin(); it != all_frag.end(); ++ it) {
    //  string ref_id; BWTIDX ref_idx;
    //  bwt.StrPosToRefPos(it->str_pos, ref_id, ref_idx);
    //  cerr << "???  " << (long long int) it->str_pos << "  " << ref_id << "  " << (long long int) ref_idx << " " << it->q_begin << " " << it->q_end << endl;
    //  cerr << "???: " << it->q_begin << " " << it->q_end << " " << it->score << endl;
    //}

    vector<FragType> valid_frag;

    int last_q_begin = -1;
    int last_q_end = -1;
    auto it_chunk_begin = all_frag.begin();
    FragType f;
    f.score = 0;
    for(auto it = all_frag.begin(); it != all_frag.end(); ++ it)
    {
        if(it->q_begin != last_q_begin || it->q_end != last_q_end)
        {
            if(f.score > 0)
            {
                for(it_chunk_begin; it_chunk_begin != it; ++ it_chunk_begin)
                {
                    // only select alignments that have the highest score
                    if(it_chunk_begin->score >= f.score)
                    {
                        TargetType t;
                        t.strand = it_chunk_begin->strand;
                        bwt.StrPosToRefPos(it_chunk_begin->str_pos, t.sid, t.begin);
                        t.end = t.begin + it_chunk_begin->str_len - 1;
                        t.num_mismatch = it_chunk_begin->num_mismatch;
                        t.cigar = new char [50];
                        strcpy(t.cigar, it_chunk_begin->cigar);
                        delete [] it_chunk_begin->cigar;
                        //cerr << "pushed solution: " << f.q_begin << " " << f.q_end << " " << t.begin << " " << t.end << " " << t.cigar << " " << it_chunk_begin->score << "  " << f.score << endl;
                        f.targets.push_back(t);
                    }
                }
                valid_frag.push_back(f);
            }
            it_chunk_begin = it;
            f.q_begin = it->q_begin;
            f.q_end = it->q_end;
            f.score = 0;
            f.targets.clear();
        }

        f.score = it->score > f.score ? it->score : f.score;
        last_q_begin = it->q_begin;
        last_q_end = it->q_end;
    }
    if(f.score > 0)
    {
        for(it_chunk_begin; it_chunk_begin != all_frag.end(); ++ it_chunk_begin)
        {
            // only select alignments that have the highest score
            if(it_chunk_begin->score >= f.score)
            {
                TargetType t;
                t.strand = it_chunk_begin->strand;
                bwt.StrPosToRefPos(it_chunk_begin->str_pos, t.sid, t.begin);
                t.end = t.begin + it_chunk_begin->str_len - 1;
                t.num_mismatch = it_chunk_begin->num_mismatch;
                t.cigar = new char [50];
                strcpy(t.cigar, it_chunk_begin->cigar);
                delete [] it_chunk_begin->cigar;
                //cerr << "pushed solution: " << f.q_begin << " " << f.q_end << " " << t.begin << " " << t.end << " " << t.cigar << " " << it_chunk_begin->score << "  " << f.score << endl;
                f.targets.push_back(t);
            }
        }
        valid_frag.push_back(f);
    }

    //cout << "valid fragment size: " << valid_frag.size() << endl;

    if(valid_frag.size() <= 0 || num_fragments <= 0)
    {
        return 0;
    }
    // using DP to find the set of best matching
    int n = valid_frag.size();
    vector<vector<int> > max_score(num_fragments);
    vector<vector<int> > prev(num_fragments);
    for(i = 0; i < num_fragments; ++ i)
    {
        max_score[i].resize(n, 0);
        prev[i].resize(n, -1);
    }
    // using DP to compute the optimal non-overlapping set for the query
    // initialization
    for(j = 0; j < n; ++ j)
    {
        //max_score[0][j] = valid_frag[j].q_end - valid_frag[j].q_begin + 1;
        max_score[0][j] = valid_frag[j].score;
    }
    // filling the rest of the matrix
    for(i = 1; i < num_fragments; ++ i)
    {
        for(j = 0; j < n; ++ j)
        {
            // this is the range of the current hit
            int score = max_score[0][j];
            // this is the lower bound score
            max_score[i][j] = max_score[i - 1][j];
            // see if the score can be higher
            for(k = 0; k < j; ++ k)
            {
                int insert_penalty = valid_frag[j].q_begin > valid_frag[k].q_end ?
                                     valid_frag[j].q_begin - valid_frag[k].q_end - 1 : valid_frag[k].q_end - valid_frag[j].q_begin + 1;
                insert_penalty = use_insert_penalty ? insert_penalty : 0;
                if(valid_frag[k].q_end < valid_frag[j].q_begin + max_overlap)
                {
                    if(score + max_score[i - 1][k] - frag_cost - insert_penalty > max_score[i][j])
                    {
                        max_score[i][j] = score + max_score[i - 1][k] - frag_cost - insert_penalty;
                        prev[i][j] = k;
                    }
                }
            }
        }
    }

    /*
    for(i = 0; i < num_fragments; ++ i) {
      for(j = 0; j < n; ++ j) {
        cerr << max_score[i][j] << "  ";
      }
      cerr << endl;
    }
    cerr << "====================" << endl;
    for(i = 0; i < num_fragments; ++ i) {
      for(j = 0; j < n; ++ j) {
        cerr << prev[i][j] << "  ";
      }
      cerr << endl;
    }
    */

    // trace-back and identify the set of sequences to be aligned
    unordered_map<int, bool> taken_pair;
    int mx = 0;
    for(j = 0; j < n; ++ j)
    {
        if(max_score[num_fragments - 1][j] > mx)
        {
            mx = max_score[num_fragments - 1][j];
        }
    }
    vector<int> current_idx;
    for(j = 0; j < n; ++ j)
    {
        if(max_score[num_fragments - 1][j] >= mx) current_idx.push_back(j);
    }

    for(k = 0; k < current_idx.size(); ++ k)
    {
        i = num_fragments - 1;
        while(current_idx[k] >= 0 && i >= 0)
        {
            FragType f;
            f = valid_frag[current_idx[k]];
            f.sol_id = k;
            selected.push_back(f);
            //cerr << "selected:  " << f.score << " " << f.q_begin << " " << f.q_end << endl;
            //for(auto it = f.targets.begin(); it != f.targets.end(); ++ it) {
            //  cerr << "selected targets:  " << it->begin << " " << it->end << " " << it->cigar << endl;
            //}
            taken_pair[current_idx[k]] = true;
            current_idx[k] = prev[i][current_idx[k]];
            -- i;
        }
    }
    // collect memory

    for(i = 0; i < valid_frag.size(); ++ i)
    {
        for(auto itf = valid_frag[i].targets.begin(); itf != valid_frag[i].targets.end(); ++ itf)
        {
            delete [] itf->cigar;
        }
    }

    //cout << "End of FindBestMatching" << endl;
    return mx;
}

