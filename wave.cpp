#include <iostream>
#include <fstream>
#include <cstdint>
#include <memory>
#include <random>
#include <chrono>
#include <bitset>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <bit>
#include <cmath>

#include "malloc_count.h"

// sdsl
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

//string filename = "resources/dblp.xml.200MB";
string filename = "resources/einstein";
//string filename = "resources/dblp.xml";
//string filename = "resources/test.txt";

uintmax_t timestamp() {
    return std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
}

/**
Wavelet Tree Version 1
*/
struct wave_bv {
    unique_ptr<bit_vector> bv;
    unique_ptr<wave_bv> left_child;
    unique_ptr<wave_bv> right_child;
    unique_ptr<rank_support_v<>> rank_bv;
    unique_ptr<vector<uint8_t>> alphabet;
    unique_ptr<bitset<256>> first_half;

    wave_bv(unique_ptr<bit_vector> bv_in,  unique_ptr<vector<uint8_t>> alphabet_in) {
        bv = move(bv_in);
        left_child = nullptr;
        right_child =  nullptr;
        rank_bv = make_unique<rank_support_v<>>(bv.get());
        alphabet = move(alphabet_in);
        first_half = make_unique<bitset<256>>();

        int sigma = (*alphabet).size();
        int split_index = sigma / 2;
        for (int i = 0; i < split_index; i++) {
            (*first_half).set((*alphabet)[i]);
        }
    }
};

int rank_wave(wave_bv& wavelet_tree, uint8_t c, int i) {
    bool is_in_first_half = (*wavelet_tree.first_half).test(c);
    if(is_in_first_half) {
        // go to left child
        int new_index = i - (*wavelet_tree.rank_bv)(i);
        if(wavelet_tree.left_child == nullptr) {
            return new_index;
        }
        return rank_wave(*wavelet_tree.left_child,  c, new_index);
    } else {
        // go to right child
        int new_index = (*wavelet_tree.rank_bv)(i);
        if(wavelet_tree.right_child == nullptr) {
            return new_index;
        }
        return rank_wave(*wavelet_tree.right_child, c, new_index);
    }
    return -1;
}

uint8_t access_wave(wave_bv& wavelet_tree, int i) {
    /*cout << "alphabet = [";
    for (uint8_t c : (*wavelet_tree.alphabet)) {
        cout << c << " ";
    }
    cout << "]" << endl;*/
    bool bit_value = (*(wavelet_tree.bv))[i];
    if(bit_value == 0) {
        // go to left child
        if(wavelet_tree.left_child == nullptr) {
            return (*wavelet_tree.alphabet)[0];
        }
        int new_index = i - (*wavelet_tree.rank_bv)(i);
        return access_wave(*wavelet_tree.left_child, new_index);
    } else {
        // go to right child
        if(wavelet_tree.right_child == nullptr) {
            return (*wavelet_tree.alphabet)[(*wavelet_tree.alphabet).size() - 1];
        }
        int new_index = (*wavelet_tree.rank_bv)(i);
        return access_wave(*wavelet_tree.right_child, new_index);
    }

    return 0;
}

unique_ptr<wave_bv> wave1(string &s) {
    // calculate alphabet size and store alphabet <-> working mapping
    int sigma = 0;
    int n = s.length();
    vector<bool> exists(255);
    vector<uint8_t> alphabet;

    for (uint8_t c : s) {
        if(!exists[c]) {
            exists[c] = true;
            sigma++;
            alphabet.push_back(c);
        }
    }

    vector<int> originalAlphabetToWorking(255);
    vector<int> workingAlphabetToOriginal(sigma);

    int workingIndex = 0;
    vector<uint8_t> alphabet_sorted;
    for (int i = 0; i < 255; i++) {
        if (exists[i]) {
            originalAlphabetToWorking[i] = workingIndex;
            workingAlphabetToOriginal[workingIndex] = i;
            workingIndex++;
            alphabet_sorted.push_back(i);
        }
    }
    //cout << "sigma = " << sigma << ", workingIndex = " << workingIndex << endl;

    // split alphabet into two bv
    unique_ptr<bit_vector> level_one = make_unique<bit_vector>(n, 0);
    string left_string;
    string right_string;
    int split_index = sigma / 2;

    for (int i = 0; i < n; i++) {
        if(originalAlphabetToWorking[(uint8_t)s[i]] >= split_index) {
            (*level_one)[i] = 1;
            right_string += s[i];
        } else {
            left_string += s[i];
        }
    }

    rank_support_v<> rank_level_one(level_one.get());
    //cout << "level_one = " << (*level_one)[0] << (*level_one)[1] << (*level_one)[2] << (*level_one)[3] << (*level_one)[4] << (*level_one)[5] << "..." << endl;
    //cout << "rank_level_one(6) = " << rank_level_one(0) << " (expected 3)" << endl;
    //cout << "size of level_one in MB: " << size_in_mega_bytes(*level_one) << endl;
    //cout << "size of rank_level_one in MB: " << size_in_mega_bytes(rank_level_one) << endl;
    //cout << "length left_string = " << left_string.length() << endl;
    //cout << "length right_string = " << right_string.length() << endl << endl;

    auto root = make_unique<wave_bv>(move(level_one), make_unique<vector<uint8_t>>(alphabet_sorted));

    // possible compute left child
    if (split_index > 1) {
        auto left_child = wave1(left_string);
        (*root).left_child = move(left_child);
    } 

    // possibly compute right child
    if (sigma > 2) {
        auto right_child = wave1(right_string);
        (*root).right_child = move(right_child);
    }

    return move(root);
}

int size_wave(wave_bv& wavelet_tree) {
    int level_size = 0;
    level_size += size_in_bytes(*wavelet_tree.bv);
    level_size += size_in_bytes(*wavelet_tree.rank_bv);
    level_size += (*wavelet_tree.alphabet).capacity() * sizeof(uint8_t);
    level_size += sizeof(bitset<256>);

    //cout << "level_size = " << level_size << endl; 

    if(wavelet_tree.left_child != nullptr) {
        level_size += size_wave(*wavelet_tree.left_child);
    }
    if(wavelet_tree.right_child != nullptr) {
        level_size += size_wave(*wavelet_tree.right_child);
    }

    return level_size;
}

/**
Efficient construction of Wavelet Trees
*/
using alphabet_char = uint8_t;
using max_occ_int = uint32_t;

void wave2(vector<alphabet_char> &T, unordered_set<alphabet_char> &alphabet) {
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

    //cout << "sigma = " << (int64_t) sigma << endl;
    //cout << "count_levels = " << (int64_t) count_levels << endl;

    /*for (alphabet_char i = 1; i < alphabet_sorted.size(); i++) {
        alphabet_char j = alphabet_sorted[i];
        histogram_alph[j] += histogram_alph[j-1];
    }*/

    /*for (alphabet_char c : alphabet_sorted) {
        cout << "histogram_alph[" << c << "] = " << histogram_alph[c] << endl;
    }*/

    // group characters into groups
    vector<max_occ_int> group_acc_frequencies;
    group_acc_frequencies.resize(sigma / 2 + 2, 0);

    vector<max_occ_int*> alph_map_to_group;
    //alph_map_to_group.reserve(biggest_char);
    cout << static_cast<size_t>(biggest_char) << endl << flush;
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
    }

    cout << 1;*/
}

void wave2(vector<alphabet_char> &T) {
    unordered_set<alphabet_char> alphabet;
    
    for (alphabet_char c : T) {
        alphabet.insert(c);
    }

    wave2(T, alphabet);
}

void sdsl_huff(string &s) {
    wt_huff<> wt;

    construct_im(wt, s, 1);

    cout << "wt[2] = " << wt[2] << endl;
    cout << "rank('a', 700) = " << wt.rank(700, 'a') << endl;
    cout << "Sigma: " << wt.sigma << endl;
}

int main(int argc, char** argv) {
    string s;

    {
        ifstream ifs(filename);
        s = string(istreambuf_iterator<char>(ifs), {});
    }
    int64_t n = s.length();
    cout << "File " << filename << " successfully loaded (n=" << n << ")" << endl << flush;

    // construction wavelet tree 1
    malloc_count_reset_peak();
    auto start_construction = timestamp();
    unique_ptr<wave_bv> wavelet_tree = wave1(s);
    auto end_construction = timestamp();
    cout << "Wavelettree for text of length " << n << " was created in " << end_construction - start_construction << " ms" << endl;

    uint64_t malloc_space_peak = malloc_count_peak();
    cout << "Wavelettree 1 malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;

    int size_wave_bytes = size_wave(*wavelet_tree);
    cout << "Größe Wavelettree = " << size_wave_bytes << " bytes (~" << size_wave_bytes / (1024 * 1024) << " MB)" << endl;


    // construction wavelet tree 2
    vector<alphabet_char> T(s.begin(), s.end());

    malloc_count_reset_peak();
    start_construction = timestamp();
    wave2(T);
    end_construction = timestamp();
    cout << "Wavelettree2 (bottom up construction) for text of length " << n << " was created in " << end_construction - start_construction << " ms" << endl;

    malloc_space_peak = malloc_count_peak();
    cout << "Wavelettree2 (bottom up construction) malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;

    cout << "Größe wave2 = " << "TODO" << " MB" << endl;

    /*uint8_t result = access_wave(*wavelet_tree, 0);
    cout << "Access[0] = " << result << " (int) " << (int) result << " (expected '<')" << endl;
    result = access_wave(*wavelet_tree, 1);
    cout << "Access[1] = " << result << " (expected '?')" << endl;
    result = access_wave(*wavelet_tree, 2);
    cout << "Access[2] = " << result << " (expected 'x')" << endl;
    result = access_wave(*wavelet_tree, 3);
    cout << "Access[3] = " << result << " (expected 'm')" << endl;

    int r = rank_wave(*wavelet_tree, '<', 60);
    cout << "wavelet_tree.rank(60, '<') = " << r << " (expected 2)" << endl;*/

    wt_huff<> wt;

    // construct sdsl huff wv tree
    malloc_count_reset_peak();
    start_construction = timestamp();
    construct_im(wt, s, 1);
    end_construction = timestamp();
    cout << "sdsl wv tree huff for text of length " << n << " was created in " << end_construction - start_construction << " ms" << endl;
    malloc_space_peak = malloc_count_peak();
    cout << "sdsl wv tree huff malloc_peak ~= " << malloc_space_peak / (1024 * 1024) << " MB" << endl;
    cout << "Größe wt_huff = " << size_in_mega_bytes(wt) << " MB" << endl;

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