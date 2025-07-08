#include <iostream>
#include <bitset>
#include <cstdint>
#include <memory>
#include <unordered_set>

// sdsl
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

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
    unique_ptr<bit_vector::select_0_type> sel_0;
    unique_ptr<bit_vector::select_1_type> sel_1;

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
        sel_0 = make_unique<bit_vector::select_0_type>(bv.get());
        sel_1 = make_unique<bit_vector::select_1_type>(bv.get());
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

unique_ptr<wave_bv> wave_old(string &s) {
    // calculate alphabet size and store alphabet <-> working mapping
    int64_t sigma = 0;
    int64_t n = s.length();
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
    int64_t split_index = sigma / 2;

    for (int64_t i = 0; i < n; i++) {
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
        auto left_child = wave_old(left_string);
        (*root).left_child = move(left_child);
    } 

    // possibly compute right child
    if (sigma > 2) {
        auto right_child = wave_old(right_string);
        (*root).right_child = move(right_child);
    }

    return move(root);
}

int size_wave(wave_bv& wavelet_tree) {
    int level_size = 0;
    level_size += size_in_bytes(*wavelet_tree.bv);
    level_size += size_in_bytes(*wavelet_tree.rank_bv);
    level_size += size_in_bytes(*wavelet_tree.sel_0);
    level_size += size_in_bytes(*wavelet_tree.sel_1);
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