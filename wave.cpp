#include <iostream>
#include <fstream>
#include <cstdint>
#include <memory>
#include <random>
#include <chrono>
#include <bitset>

// sdsl
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

string filename = "resources/dblp.xml";

uintmax_t timestamp() {
    return std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
}

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

    // construction wavelet tree
    auto start_construction = timestamp();
    unique_ptr<wave_bv> wavelet_tree = wave1(s);
    auto end_construction = timestamp();
    cout << "Wavelettree for text of length " << n << " was created in " << end_construction - start_construction << " ms" << endl;

    uint8_t result = access_wave(*wavelet_tree, 0);
    cout << "Access[0] = " << result << " (int) " << (int) result << " (expected '<')" << endl;
    result = access_wave(*wavelet_tree, 1);
    cout << "Access[1] = " << result << " (expected '?')" << endl;
    result = access_wave(*wavelet_tree, 2);
    cout << "Access[2] = " << result << " (expected 'x')" << endl;
    result = access_wave(*wavelet_tree, 3);
    cout << "Access[3] = " << result << " (expected 'm')" << endl;

    int r = rank_wave(*wavelet_tree, '<', 60);
    cout << "wavelet_tree.rank(60, '<') = " << r << " (expected 2)" << endl;

    wt_huff<> wt;

    // construct sdsl huff wv tree
    start_construction = timestamp();
    construct_im(wt, s, 1);
    end_construction = timestamp();
    cout << "sdsl wv tree huff for text of length " << n << " was created in " << end_construction - start_construction << " ms" << endl;

    // create test indices
    const int m = 1000000;
    vector<int> indices;

    random_device rd; 
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, n - 1);
    for (int64_t i = 0; i < m; i++) {
        indices.push_back(dist(gen));
    }

    // assert correct
    for (int i = 0; i < indices.size() && i < 5; i++) {
        int j = indices[i];
        cout << "wavelet_tree[" << j << "] = " << access_wave(*wavelet_tree, j)
            << ", wt_huff[" << j << "] = " << wt[j] << endl;
    }

    uint8_t rank_char = 'x';
    for (int i = 0; i < indices.size() && i < 5; i++) {
        int j = indices[i];
        cout << "wavelet_tree.rank(" << j << ", '" << rank_char << "') = " << rank_wave(*wavelet_tree, rank_char, j)
            << ", wt_huff.rank(" << j << ", '" << rank_char << "') = " << wt.rank(j, rank_char) << endl;
    }

    // test speed access
    auto start = timestamp();
    int64_t x = 0;
    for (int i : indices) {
        x += (int)wt[i];
    }
    auto end = timestamp();
    cout << "SDSL WT Huff Access with " << m << " iterations finished in " << end - start << " ms" << endl;

    start = timestamp();
    x = 0;
    for (int i : indices) {
        x += access_wave(*wavelet_tree, i);
    }
    end = timestamp();
    cout << "Wavelettree Access with " << m << " iterations finished in " << end - start << " ms" << endl;
}