#include <iostream>
#include <fstream>
#include <cstdint>
#include <memory>

// sdsl
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

string filename = "resources/dblp.xml";



struct wave_bv {
    unique_ptr<bit_vector> bv;
    unique_ptr<wave_bv> left_child;
    unique_ptr<wave_bv> right_child;
    unique_ptr<rank_support_v<>> rank_bv;
    unique_ptr<vector<uint8_t>> alphabet;

    wave_bv(unique_ptr<bit_vector> bv_in,  unique_ptr<vector<uint8_t>> alphabet_in) {
        bv = move(bv_in);
        left_child = nullptr;
        right_child =  nullptr;
        rank_bv = make_unique<rank_support_v<>>(bv.get());
        alphabet = move(alphabet_in);
    }
};

uint8_t access_wave(wave_bv& wavelet_tree, int i) {
    cout << "alphabet = [";
    for (uint8_t c : (*wavelet_tree.alphabet)) {
        cout << c << " ";
    }
    cout << "]" << endl;

    cerr << i << endl;
    if(wavelet_tree.bv == nullptr) {
            cerr << "INVALID NULL POINTER" << endl << flush;
            return 0;
        }
    bool bit_value = (*(wavelet_tree.bv))[i];
    cout << "bit_value[" << i << "] = " << bit_value << endl;
    if(bit_value == 0) {
        // go to left child
        if(wavelet_tree.left_child == nullptr) {
            return (*wavelet_tree.alphabet)[0];
        }
        if(wavelet_tree.rank_bv == nullptr) {
            cerr << "INVALID NULL POINTER" << endl << flush;
            return 0;
        } else {
            //rank_support_v<> rank_temp(wavelet_tree.bv.get());
            int new_index = i - (*wavelet_tree.rank_bv)(i);
            cout << "new_index = " << new_index << endl;
            if(wavelet_tree.left_child == nullptr) {
                cerr << "INVALID NULL POINTER" << endl << flush;
                return 0;
            }
            return access_wave(*wavelet_tree.left_child, new_index);
        }
    } else {
        // go to right child
        if(wavelet_tree.right_child == nullptr) {
            cout << "hello world" << endl;
            return (*wavelet_tree.alphabet)[(*wavelet_tree.alphabet).size() - 1];
        }
        int new_index = (*wavelet_tree.rank_bv)(i);
        cout << "new_index_right = " << new_index << endl;
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
    cout << "sigma = " << sigma << ", workingIndex = " << workingIndex << endl;

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
    cout << "level_one = " << (*level_one)[0] << (*level_one)[1] << (*level_one)[2] << (*level_one)[3] << (*level_one)[4] << (*level_one)[5] << "..." << endl;
    cout << "rank_level_one(6) = " << rank_level_one(0) << " (expected 3)" << endl;
    cout << "size of level_one in MB: " << size_in_mega_bytes(*level_one) << endl;
    cout << "size of rank_level_one in MB: " << size_in_mega_bytes(rank_level_one) << endl;
    cout << "length left_string = " << left_string.length() << endl;
    cout << "length right_string = " << right_string.length() << endl << endl;

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

void sdsl_huff(string &s) {
    wt_huff<> wt;

    construct_im(wt, s, 1);

    cout << "wt[2] = " << wt[2] << endl;
    cout << "rank('a', 700) = " << wt.rank(700, 'a') << endl;
    cout << "Sigma: " << wt.sigma << endl;
}

int main(int argc, char** argv) {
    cout << "Hi, was geht" << endl;

    string s;

    {
        ifstream ifs(filename);
        s = string(istreambuf_iterator<char>(ifs), {});
    }
    int64_t n = s.length();
    cout << "File " << filename << " successfully loaded (n=" << n << ")" << endl << flush;

    unique_ptr<wave_bv> wavelet_tree = wave1(s);

    uint8_t result = access_wave(*wavelet_tree, 0);
    cout << "Access[0] = " << result << " (int) " << (int) result << " (expected '<')" << endl;
    result = access_wave(*wavelet_tree, 1);
    cout << "Access[1] = " << result << " (expected '?')" << endl;
    result = access_wave(*wavelet_tree, 2);
    cout << "Access[2] = " << result << " (expected 'x')" << endl;
    result = access_wave(*wavelet_tree, 3);
    cout << "Access[3] = " << result << " (expected 'm')" << endl;

    sdsl_huff(s);
}