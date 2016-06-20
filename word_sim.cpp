#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#define random(x) (rand()%x)
#define ONE 0.999

using namespace std;

vector<pair<int, double> > get_repr(string s) {
    vector<double> mid_ret;
    int left = 0, right = 0;
    for (; right < s.length(); right++) {
        if (s[right] == '<') return vector<pair<int, double> >();
        if (s[right] == '\t') {
            double attr = atof(s.substr(left, right - left).c_str());
            left = right+1;
            right = left-1;
            mid_ret.push_back(attr);
        }
    }
    double attr = atof(s.substr(left, right - left).c_str());
    mid_ret.push_back(attr);
    
    vector<pair<int, double> > ret;
    for (int i = 1; i < mid_ret.size(); i += 2) {
        ret.push_back(make_pair(int(mid_ret[i]), mid_ret[i+1]));
    }
    sort(ret.begin(), ret.end());
    return ret;
}

double cosine_sim(const vector<pair<int, double> > &repr1, const vector<pair<int, double> > &repr2) {
    if (repr1.size() != repr2.size()) return 0;
    double m1 = 0.0, m2 = 0.0, dot_product = 0.0;
    int i = 0, j = 0;
    for (; i < repr1.size() && j < repr2.size();) {
        if (repr1[i].first < repr2[j].first) {
            m1 += repr1[i].second * repr1[i].second;
            i++;
        }
        else if (repr1[i].first > repr2[j].first) {
            m2 += repr2[j].second * repr2[j].second;
            j++;
        }
        else {
            dot_product += repr1[i].second * repr2[j].second;
            m1 += repr1[i].second * repr1[i].second;
            m2 += repr2[j].second * repr2[j].second;
            i++;j++;
        }
    }
    while (i < repr1.size()) {
        m1 += repr1[i].second * repr1[i].second;
        i++;
    }
    while (j < repr2.size()) {
        m2 += repr2[j].second * repr2[j].second;
        j++;
    }
    m1 = sqrt(m1); m2 = sqrt(m2);
    return dot_product / (m1 * m2);
}

int main(int argc, char *argv[]) {
    string dict_file_name = argv[1]; //"clusters300k_old.txt";
    string word_file_name = argv[2];
    string output_dir = argv[3];
    //string output_file_name = "/Users/danielwu/Documents/coding/singluariti/base/data/semi_markov/sim_data/clusters300k_sim.txt";
    ifstream dict_in_file(dict_file_name.c_str());
    ifstream word_in_file(word_file_name.c_str());
    //ofstream out_file(output_file_name.c_str());
    string line;
    unordered_map<string, vector<pair<int, double> > > mp;
    unordered_map<string, vector<pair<int, double> > > dict_mp;
    vector<string> voca;
    vector<string> dict;
    unordered_set<string> sw_dict;
    
    //getline(in_file, line);
    while (getline(dict_in_file,line)) {
        int idx = 0;
        for (; idx < line.length(); idx++) {
            if (line[idx] == '-')
                break;
        }
        //cout << idx << endl;
        string word = line.substr(0, idx);
        // cout << word << " " << word.length() << endl;
        string word_sim;
        int step = 3;
        for (int i = 0; i < word.length();) {
            if (step == 3)
                word_sim += word.substr(i,step);
            i += step;
            if (step == 3) step = 1;
            else step = 3;
        }
        //cout << word_sim << word_sim.length() << endl;
        
        vector<pair<int, double> > repr = get_repr(line.substr(idx));
        if (repr.size() >= 8) {
            dict.push_back(word_sim);
            //for (int i = 0; i < repr.size(); i++)
            //cout << repr[i] << endl;
            
            dict_mp[word_sim] = repr;
        }
        if (repr.size() == 0) {
            sw_dict.insert(word_sim);
        }
    }

   while (getline(word_in_file,line)) {
        int idx = 0;
        for (; idx < line.length(); idx++) {
            if (line[idx] == '-')
                break;
        }
        //cout << idx << endl;
        string word = line.substr(0, idx);
        // cout << word << " " << word.length() << endl;
        string word_sim;
        int step = 3;
        for (int i = 0; i < word.length();) {
            if (step == 3)
                word_sim += word.substr(i,step);
            i += step;
            if (step == 3) step = 1;
            else step = 3;
        }
        //cout << word_sim << word_sim.length() << endl;
        
        vector<pair<int, double> > repr = get_repr(line.substr(idx));
        if (repr.size() >= 8) {
            voca.push_back(word_sim);
            //for (int i = 0; i < repr.size(); i++)
            //cout << repr[i] << endl;
            
            mp[word_sim] = repr;
        }
    }

    unordered_set<int> list;
    for (int i = 0; i < 10000; i++) {
        list.insert(random(290000));
    }
    int cnt = 0; 
    for (int i = 0; i < voca.size(); i++) {
        //int cnt = 0;
//        if (list.find(i) == list.end()) continue;
        //if (voca[i] != "共同的目标") continue;
        //if (i != 100) continue;
	if (i >= 20000) continue;
        if (voca[i].length() <= 6) continue;
        string out_file_name = output_dir + "/" + word_file_name.substr(0, word_file_name.length() - 4) + "_" + voca[i] + ".txt";
        ofstream out_file(out_file_name.c_str());
        //out_file << voca[i] << " ";
        //bool is_print = false;
        // output the similarity with its prefix and suffix
	for (int k = 3; k < voca[i].length(); k += 3) {
	    string s = voca[i].substr(0,k);
	    if (dict_mp.find(s) != dict_mp.end()) {
	        out_file << voca[i] << "\t" << s << "\t" << cosine_sim(mp[voca[i]],dict_mp[s]) << endl;
	    }
	}

	for (int k = voca[i].length() - 3; k > 0; k -= 3) {
	    string s = voca[i].substr(k);
	    if (dict_mp.find(s) != dict_mp.end()) {
	        out_file << voca[i] << "\t" << s << "\t" << cosine_sim(mp[voca[i]], dict_mp[s]) << endl;
	    }   
	}
	out_file << endl;

        vector<pair<double, int> > indexs;
        vector<vector<string> > phrase_pairs;
        for (int j = 0; j < dict.size(); j++) {
            //cout << voca[j] << voca[j].length() << endl;
            if (voca[i] == dict[j] || dict[j].length() <= 6) continue;
            //if (voca[j] != "美好的人生") continue;
            double tot_sim = cosine_sim(mp[voca[i]], dict_mp[dict[j]]);
            if (tot_sim <= 0.2) continue;
            //if (voca[j].length() == 3) continue;
            //out_file << voca[i] << " " << voca[j] << " ";
            //           out_file << voca[j] << " ";
            string voca_i = voca[i];
            string voca_j = dict[j];
            vector<vector<double> > dp(30, vector<double>(30, 0.001));
            vector<vector<pair<int, int>> > pos(30, vector<pair<int, int>>(30,{0,0}));
            for (int k = 2; k < voca_i.size(); k += 3) {
                for (int l = 2; l < voca_j.size(); l += 3) {
                    if (k == voca_i.size()-1 && l != voca_i.size()-1) continue;
                    if (k != voca_i.size()-1 && l == voca_j.size()-1) continue;
                    pos[k][l] = make_pair(k, l);
                    if (voca_i.substr(0,k+1) == voca_j.substr(0,l+1) && dict_mp.find(voca_i.substr(0,k+1)) != dict_mp.end()) {
                        if (k!= voca_i.size()-1 && l != voca_j.size()-1) {
                            dp[k][l] = 1;
                        }
                        continue;
                    }
                    if (dict_mp.find(voca_i.substr(0,k+1)) != dict_mp.end() && dict_mp.find(voca_j.substr(0,l+1)) != dict_mp.end()) {
                        string str1 = voca_i.substr(0,k+1);
                        string str2 = voca_j.substr(0,l+1);
                        vector<pair<int, double> > repr1 = dict_mp[str1];
                        vector<pair<int, double> > repr2 = dict_mp[str2];
                        
                        if (k != voca_i.size()-1 && l != voca_j.size()-1) {
                            double sim = cosine_sim(repr1, repr2);
                            dp[k][l] = sim;
                            //if (sim > 0.9999)
                            //    continue;
                        }
                    }
                    for (int str1 = k-2; str1 > 0; str1 -= 3) {
                        for (int str2 = l-2; str2 > 0; str2 -= 3) {
                            string tmp_str1 = voca_i.substr(str1, k - str1 + 1);
                            string tmp_str2 = voca_j.substr(str2, l - str2 + 1);
                            string tmp_str3 = voca_i.substr(0, str1);
                            string tmp_str4 = voca_j.substr(0, str2);
                            
                            if (sw_dict.find(tmp_str1) != sw_dict.end() && sw_dict.find(tmp_str2) != sw_dict.end()) {
                                if (dp[k][l] < dp[str1-1][str2-1]*ONE) {
                                    dp[k][l] = dp[str1-1][str2-1]*ONE;
                                    pos[k][l] = make_pair(str1-1, str2-1);
                                }
                            }
                            else {
                                if (tmp_str1 == tmp_str2) {
                                    if (dp[str1-1][str2-1]*ONE > dp[k][l]) {
                                        dp[k][l] = dp[str1-1][str2-1]*ONE;
                                        pos[k][l] = make_pair(str1-1, str2-1);
                                        continue;
                                    }
                                }
                            
                                if (dict_mp.find(tmp_str1) != dict_mp.end() && dict_mp.find(tmp_str2) != dict_mp.end()) {
                                    if (dp[str1-1][str2-1] * cosine_sim(dict_mp[tmp_str1], dict_mp[tmp_str2]) > dp[k][l]) {
                                        dp[k][l] = dp[str1-1][str2-1] * cosine_sim(dict_mp[tmp_str1], dict_mp[tmp_str2]);
                                        pos[k][l] = make_pair(str1-1, str2-1);
                                    }
                                }
                            }
                        }
                    }
                    
                }
            }
            
            string phrase1 = voca_i;
            string phrase2 = voca_j;
            
            int phrase1_ptr = voca_i.length()-1;
            int phrase2_ptr = voca_j.length()-1;
            
            int cnt = 0;
            
            do{
                phrase1.insert(phrase1_ptr+1, "|");
                phrase2.insert(phrase2_ptr+1, "|");
                //cnt++;
                
                int tmp_phrase1_ptr = pos[phrase1_ptr][phrase2_ptr].first;
                int tmp_phrase2_ptr = pos[phrase1_ptr][phrase2_ptr].second;
                
                if (tmp_phrase2_ptr != phrase2_ptr || tmp_phrase1_ptr != phrase1_ptr)
                    cnt++;
                
                phrase1_ptr = tmp_phrase1_ptr;
                phrase2_ptr = tmp_phrase2_ptr;
            }
            while (phrase1_ptr != pos[phrase1_ptr][phrase2_ptr].first && phrase2_ptr != pos[phrase1_ptr][phrase2_ptr].second);
            
            phrase1.insert(phrase1_ptr+1, "|");
            phrase2.insert(phrase2_ptr+1, "|");
            
            // if (cnt > 0) {
            
            //cout << dp[voca[i].size()-1][voca[j].size()-1] << endl;
            //if (dp[voca[i].size()-1][voca[j].size()-1] > 0.5) cnt++;
            vector<string> phrase_pair;
            phrase_pair.push_back(phrase1);
            phrase_pair.push_back(phrase2);
            phrase_pairs.push_back(phrase_pair);
            
            pair<double, int> index{dp[voca_i.size()-1][voca_j.size()-1]*tot_sim, phrase_pairs.size()-1};
            indexs.push_back(index);
            
            //is_print = true;
            //   }
            //}
        }
        sort(indexs.begin(), indexs.end());
        for (int i = indexs.size() - 1; i >= 0; i--) {
            if (indexs[i].first < 0.1) break;
	    out_file << phrase_pairs[indexs[i].second][0] << " " << phrase_pairs[indexs[i].second][1] << " " << indexs[i].first << endl;
        }
        // if (is_print)
        //   out_file << endl;
    }
   // cout << "cnt = " << cnt << endl;
    // cout << longest_phrase_idx << " " << longest_phrase_length << endl;
    
    // cout << longest_phrase_idx << " " << longest_phrase_length << endl;
}
