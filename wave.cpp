#include <iostream>
#include <fstream>
#include <cstdint>
#include <memory>
#include <random>
#include <chrono>
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>

#include "malloc_count.h"

#include "wave-old.hpp"

// sdsl
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

string filename = "resources/english";
//string filename = "resources/dblp.xml.200MB";
//string filename = "resources/einstein";
//string filename = "resources/dblp.xml";
//string filename = "resources/test.txt";

uintmax_t timestamp() {
    return std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
}

/**
Bottom-up construction of Wavelet Trees
*/
using alphabet_char = uint8_t;
using max_occ_int = uint32_t;

struct wv_tree {
    vector<unique_ptr<bit_vector>> levels;
    vector<unique_ptr<rank_support_v<>>> rank_levels;
    vector<unique_ptr<bit_vector::select_0_type>> sel_0_levels;
    vector<unique_ptr<bit_vector::select_1_type>> sel_1_levels;
};

wv_tree wave(vector<alphabet_char> &T, unordered_set<alphabet_char> &alphabet) {
    auto preprocessing_start = timestamp();

    alphabet_char sigma = alphabet.size();
    alphabet_char count_levels = bit_width(static_cast<uint64_t>(sigma - 1));
    max_occ_int n = T.size();

    // calculate biggest char of alphabet
    alphabet_char biggest_char = 0;
    for (alphabet_char c : alphabet) {
        if (c > biggest_char) {
            biggest_char = c;
        }
    }

    vector<max_occ_int> histogram_alph;
    histogram_alph.resize(static_cast<size_t>(biggest_char + 1), 0);

    // update histogram map to frequencies of the symbols in the actual text T
    for (alphabet_char c : T) {
        histogram_alph[c] += 1;
    }

    vector<alphabet_char> alphabet_sorted(alphabet.begin(), alphabet.end());
    sort(alphabet_sorted.begin(), alphabet_sorted.end());

    auto other_prep_start = timestamp();

    vector<alphabet_char> alph_map_to_sorted_index;
    alph_map_to_sorted_index.resize(static_cast<size_t>(biggest_char + 1), 0);

    for (alphabet_char i = 0; i < sigma; i++) {
        alph_map_to_sorted_index[alphabet_sorted[i]] = i;
    }

    // group characters into groups
    vector<max_occ_int> group_acc_frequencies;
    group_acc_frequencies.resize(sigma / 2 + 2, 0);

    vector<max_occ_int*> alph_map_to_group;
    alph_map_to_group.resize(static_cast<size_t>(biggest_char + 1), nullptr);


    // building the levels
    vector<unique_ptr<bit_vector>> levels;
    levels.reserve(count_levels);
    for (alphabet_char i = 0; i < count_levels; i++) {
        levels.push_back(make_unique<bit_vector>(n, 0));
    }

    cout << "duration preprocessing " << timestamp() - preprocessing_start << " ms" << endl << flush;

    auto start_main = timestamp();
    alphabet_char level_index = count_levels;
    do {
        level_index--;
        alphabet_char reverse_level_index = count_levels - level_index;
        uint64_t group_size = pow(2, reverse_level_index);
        cout << "group_size = " << (int) group_size << endl;

        alphabet_char group_index = 0;
        group_acc_frequencies[0] = 0;
        for (alphabet_char i = 0; i < alphabet_sorted.size(); i++) {
            alphabet_char c = alphabet_sorted[i];
            if (i % group_size == 0 && i > 0) {
                group_index++;
            }
            if (i % group_size == 0 && i / group_size < sigma / 2 + 1) {
                group_acc_frequencies[i / group_size + 1] = group_acc_frequencies[i / group_size];
            }
            if (i / group_size < sigma / 2 + 1) {
                group_acc_frequencies[i / group_size + 1] += histogram_alph[c];
            }
            alph_map_to_group[c] = &(group_acc_frequencies[group_index]);
        }

        /*for (alphabet_char c : alphabet_sorted) {
            cout << "*alph_map_to_group[" << c << "] = " << *alph_map_to_group[c] << endl;
        }*/

        auto& level = *levels[level_index];

        vector<bool> bv_values;
        bv_values.resize(biggest_char + 1, 0);
        for (alphabet_char c : alphabet) {
            bv_values[c] = (alph_map_to_sorted_index[c] & (group_size - 1)) >= (group_size >> 1); // äquivalent zu: alph_map_to_sorted_index[c] % group_size >= group_size / 2
        }

        auto inner_loop = timestamp();

        for (max_occ_int i = 0; i < n; i++) {
            alphabet_char c = T[i];
            max_occ_int bv_index = (*alph_map_to_group[c])++;
            level[bv_index] = bv_values[c];
        }

        cout << "duration inner loop " << (int) level_index << ": " << timestamp() - inner_loop << " ms" << endl;
    } while (level_index > 0);

    auto end_main = timestamp();
    cout << "time main loop " << end_main - start_main << " ms" << endl;

    auto start_rank = timestamp();

    vector<unique_ptr<rank_support_v<>>> rank_levels;
    for (alphabet_char level_index = 0; level_index < count_levels; level_index++) {
        rank_levels.push_back(make_unique<rank_support_v<>>(levels[level_index].get()));
    }

    cout << "duration build rank support " << timestamp() - start_rank << " ms" << endl;

    auto start_select = timestamp();

    vector<unique_ptr<bit_vector::select_0_type>> sel_0_levels;
    vector<unique_ptr<bit_vector::select_1_type>> sel_1_levels;

    for (alphabet_char level_index = 0; level_index < count_levels; level_index++) {
        sel_0_levels.push_back(make_unique<bit_vector::select_0_type>(levels[level_index].get()));
        sel_1_levels.push_back(make_unique<bit_vector::select_1_type>(levels[level_index].get()));
    }

    cout << "duration build select support " << timestamp() - start_select << " ms" << endl;

    // print wavelettree
    /*for (int i = 0; i < count_levels; i++) {
        for (int j = 0; j < n; j++) {
            cout << (*levels[i])[j];
        }
        cout << endl;
    }*/

    return {move(levels), move(rank_levels), move(sel_0_levels), move(sel_1_levels)};
}

wv_tree wave(vector<alphabet_char> &T) {
    unordered_set<alphabet_char> alphabet;
    
    for (alphabet_char c : T) {
        alphabet.insert(c);
    }

    return wave(T, alphabet);
}

size_t size_wave(wv_tree &wt) {
    size_t wt_size = 0;
    for (int i = 0; i < wt.levels.size(); i++) {
        wt_size += size_in_bytes(*wt.levels[i]);
        wt_size += size_in_bytes(*wt.rank_levels[i]);
        wt_size += size_in_bytes(*wt.sel_0_levels[i]);
        wt_size += size_in_bytes(*wt.sel_1_levels[i]);
    }
    return wt_size;
}

int main(int argc, char** argv) {
    string s;

    {
        ifstream ifs(filename);
        s = string(istreambuf_iterator<char>(ifs), {});
    }
    int64_t n = s.length();
    cout << "File " << filename << " successfully loaded (n=" << n << ")" << endl << endl << flush;


    // bottom-up construction wavelet tree
    cout << "===bottom-up construction of wavelet tree===" << endl;
    vector<alphabet_char> T(s.begin(), s.end());

    malloc_count_reset_peak();
    auto start_construction = timestamp();
    wv_tree wt2 = wave(T);
    auto end_construction = timestamp();
    cout << "-> bottom-up construction of wavelet tree (with n=" << n << ") was created in " << end_construction - start_construction << " ms" << endl;

    uint64_t malloc_space_peak = malloc_count_peak();
    cout << "-> malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;
    size_t wt_size = size_wave(wt2);
    cout << "-> size of wavelet tree = " << wt_size << " bytes (~ " << wt_size / (1024 * 1024) << " MB)" << endl << endl;

    // old recursive construction wavelet tree
    cout << "===old recursive construction of wavelet tree===" << endl;
    malloc_count_reset_peak();
    start_construction = timestamp();
    unique_ptr<wave_bv> wavelet_tree = wave_old(s);
    end_construction = timestamp();
    cout << "-> old recursive construction of wavelet tree (with n=" << n << ") was created in " << end_construction - start_construction << " ms" << endl;

    malloc_space_peak = malloc_count_peak();
    cout << "-> malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;

    int size_wave_bytes = size_wave(*wavelet_tree);
    cout << "-> size of wavelet tree = " << size_wave_bytes << " bytes (~" << size_wave_bytes / (1024 * 1024) << " MB)" << endl << endl;

    // construct sdsl huff wv tree
    cout << "===sdsl construction of huffman-shaped wavelet tree===" << endl;
    wt_huff<> wt;
    malloc_count_reset_peak();
    start_construction = timestamp();
    construct_im(wt, s, 1);
    end_construction = timestamp();
    cout << "-> sdsl construction of huffman-shaped wavelet tree (with n=" << n << ") was created in " << end_construction - start_construction << " ms" << endl;
    malloc_space_peak = malloc_count_peak();
    cout << "-> malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;
    cout << "-> size of wavelet tree = " << size_in_mega_bytes(wt) << " MB" << endl << endl;

    // construct sdsl int balanced wv tree
    /*wm_int<> wm_int_bal;

    malloc_count_reset_peak();
    start_construction = timestamp();
    construct_im(wm_int_bal, s, 1);
    end_construction = timestamp();
    cout << "sdsl wm tree integer balanced for text of length " << n << " was created in " << end_construction - start_construction << " ms" << endl;
    malloc_space_peak = malloc_count_peak();
    cout << "sdsl wm tree integer balanced malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;
    cout << "Größe wm_int_bal = " << size_in_mega_bytes(wm_int_bal) << " MB" << endl << endl;*/

    // construct sdsl byte balanced wv tree
    cout << "===sdsl construction of balanced wavelet tree===" << endl;
    wt_blcd<> wt_byte_bal;

    malloc_count_reset_peak();
    start_construction = timestamp();
    construct_im(wt_byte_bal, s, 1);
    end_construction = timestamp();
    cout << "-> sdsl construction of balanced wavelet tree (with n=" << n << ") was created in " << end_construction - start_construction << " ms" << endl;
    malloc_space_peak = malloc_count_peak();
    cout << "-> malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;
    cout << "-> size of wavelet tree = " << size_in_mega_bytes(wt_byte_bal) << " MB" << endl << endl;

    // construct sdsl hu-tucker wv tree
    cout << "===sdsl construction of hu-tucker shaped wavelet tree===" << endl;
    wt_hutu<> wt_hutucker;

    malloc_count_reset_peak();
    start_construction = timestamp();
    construct_im(wt_hutucker, s, 1);
    end_construction = timestamp();
    cout << "-> sdsl construction of hu-tucker shaped wavelet tree (with n=" << n << ") was created in " << end_construction - start_construction << " ms" << endl;
    malloc_space_peak = malloc_count_peak();
    cout << "-> malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;
    cout << "-> size of wavelet tree = " << size_in_mega_bytes(wt_hutucker) << " MB" << endl << endl;

    // construct sdsl int balanced wv tree
    /*wt_int<> wt_int_bal;

    malloc_count_reset_peak();
    start_construction = timestamp();
    construct_im(wt_int_bal, s, 1);
    end_construction = timestamp();
    cout << "sdsl wv tree integer balanced for text of length " << n << " was created in " << end_construction - start_construction << " ms" << endl;
    malloc_space_peak = malloc_count_peak();
    cout << "sdsl wv tree integer balanced malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;
    cout << "Größe wt_int_bal = " << size_in_mega_bytes(wt_int_bal) << " MB" << endl << endl;*/


    // create test indices
    const int m = 1000000;
    vector<int> indices;

    random_device rd; 
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, n - 1);
    for (int64_t i = 0; i < m; i++) {
        indices.push_back(dist(gen));
    }

    int sigma = 0;
    vector<bool> exists(255);
    vector<uint8_t> alphabet;

    for (uint8_t c : s) {
        if(!exists[c]) {
            exists[c] = true;
            sigma++;
            alphabet.push_back(c);
        }
    }
    vector<uint8_t> alphabet_sorted;
    for (int i = 0; i < 255; i++) {
        if (exists[i]) {
            alphabet_sorted.push_back(i);
        }
    }

    vector<uint8_t> random_chars;
    uniform_int_distribution<int>dist_alph(0, sigma - 1);
    for (int64_t i = 0; i < m; i++){
        random_chars.push_back(dist_alph(gen));
    }
    

    // assert correct
    for (int i = 0; i < indices.size() && i < 5; i++) {
        int j = indices[i];
        cout << "wavelet_tree[" << j << "] = " << access_wave(*wavelet_tree, j)
            << ", wt_huff[" << j << "] = " << wt[j] << endl;
    }

    for (int i = 0; i < indices.size() && i < 5; i++) {
        int j = indices[i];
        uint8_t rank_char = random_chars[i];
        cout << "wavelet_tree.rank(" << j << ", '" << rank_char << "') = " << rank_wave(*wavelet_tree, rank_char, j)
            << ", wt_huff.rank(" << j << ", '" << rank_char << "') = " << wt.rank(j, rank_char) << endl;
    }

    // test speed access
    auto start = timestamp();
    int64_t x = 0;
    for (int i : indices) {
        x += access_wave(*wavelet_tree, i);
    }
    auto end = timestamp();
    cout << "Wavelettree Access with " << m << " iterations finished in " << end - start << " ms" << endl;

    start = timestamp();
    x = 0;
    for (int i : indices) {
        x += (int)wt[i];
    }
    end = timestamp();
    cout << "SDSL WT Huff Access with " << m << " iterations finished in " << end - start << " ms" << endl;

    // test speed rank
    start = timestamp();
    x = 0;
    for (int i = 0; i < m; i++) {
        int j = indices[i];
        uint8_t rank_char = random_chars[i];
        x += rank_wave(*wavelet_tree, rank_char, j);
    }
    end = timestamp();
    cout << "Wavelettree Rank with " << m << " iterations finished in " << end - start << " ms" << endl;

    start = timestamp();
    x = 0;
    for (int i = 0; i < m; i++) {
        int j = indices[i];
        uint8_t rank_char = random_chars[i];
        x += wt.rank(j, rank_char);
    }
    end = timestamp();
    cout << "SDSL WT Huff Rank with " << m << " iterations finished in " << end - start << " ms" << endl;
}