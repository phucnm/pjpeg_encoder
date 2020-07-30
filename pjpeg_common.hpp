//
//  pjpeg_common.cpp
//  Toy Image Compressor
//
//  Created by Phuc Nguyen on 07/21/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#ifndef PJPEG_COMMON_HPP
#define PJPEG_COMMON_HPP

#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

const int LOW_QUALITY_CODE = 0;
const int MEDIUM_QUALITY_CODE = 1;
const int HIGH_QUALITY_CODE = 2;

const vector<vector<int>> const_medium_quantize_vector{
    {
        {16,  11,  10,  16,  24,  40,  51,  61},
        {12,  12,  14,  19,  26,  58,  60,  55},
        {14,  13,  16,  24,  40,  57,  69,  56},
        {14,  17,  22,  29,  51,  87,  80,  62},
        {18,  22,  37,  56,  68, 109, 103,  77},
        {24,  35,  55,  64,  81, 104, 113,  92},
        {49,  64,  78,  87, 103, 121, 120, 101},
        {72,  92,  95,  98, 112, 100, 103,  99}
    }
};

const vector<vector<int>> zigzag {
    {
        { 0 , 1 ,  5,  6, 14, 15, 27, 28},
        { 2 , 4 ,  7, 13, 16, 26, 29, 42},
        { 3 , 8 , 12, 17, 25, 30, 41, 43},
        { 9 , 11, 18, 24, 31, 40, 44, 53},
        { 10, 19, 23, 32, 39, 45, 52, 54},
        { 20, 22, 33, 38, 46, 51, 55, 60},
        { 21, 34, 37, 47, 50, 56, 59, 61},
        { 35, 36, 48, 49, 57, 58, 62, 63},
    }
};

// A hash function used to hash a pair of any kind
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

inline int next_mul(int n, int mod) {
    if (n % mod == 0) {
        return n;
    }
    return n + (mod - (n % mod));
}

/*
 Helper method to reverse bit in the same num_bits range
 */
inline u32 reverse_bits(u32 val, int num_bits) {
    u32 res = 0;
    u32 _val = val;
    int remain = num_bits;
    while (_val) {
        res <<= 1;
        res |= _val & 1;
        _val >>= 1;
        remain--;
    }
    res <<= remain;
    return res;
}

//Generate a map of symbol: pair <codelen, code word>
//Symbol can be negative due to DCT

template<typename T>
unordered_map<T, pair<u16, u16>> code_lengths_to_code_table(const unordered_map<T, u32>& lengths) {
    u32 max_len = 0;
    for (auto l: lengths) {
        if (l.second > max_len)
            max_len = l.second;
    }
    vector<u32> bl_count(max_len+1, 0);
    //Step 1
    for (auto l: lengths) {
        bl_count[l.second]++;
    }

    //Step 2
    u32 code = 0;
    bl_count[0] = 0;
    vector<u32> next_code(max_len+1, 0);
    for (u32 bits = 1; bits <= max_len; bits++) {
        code = (code + bl_count[bits-1]) << 1;
        next_code[bits] = code;
    }

    //Step 3
    //Sort the lengths to make sure reproducible
    vector<T> syms;
    syms.reserve(lengths.size());
    for (auto l: lengths)
        syms.push_back(l.first);
    sort(syms.begin(), syms.end());
    unordered_map<T, pair<u16, u16>> code_table;
    for (T sym: syms) {
        int len = lengths.at(sym);
        code_table[sym] = make_pair(len, next_code[len]);
        next_code[len]++;
    }
    return code_table;
}

inline int bits_count(int n) {
    if (n==0)
        return 0;
    if (n < 0)
        n = -n;
    return floor(log2((double)n)) + 1;
}

struct RunSize {
    int run, size;
    bool operator==(const RunSize &other) const {
        return run == other.run
        && size == other.size;
    }
    bool operator<(const RunSize& rhs) const {
        if (run == rhs.run)
            return size < rhs.size;
        return run < rhs.run;

    }
};

namespace std {
    template<> struct hash<RunSize> {
        std::size_t operator()(const RunSize& rs) const {
            return ((hash<int>()(rs.run)
                     ^ (hash<int>()(rs.size) << 1)));
        }
    };
}

inline vector<pair<RunSize, int>> rle_encoding(const vector<int>& input) {
    vector<pair<RunSize, int>> encoded;
    int zero_count = 0;
    for (unsigned int i = 0; i < input.size(); i++) {
        if (input[i] != 0) {
            encoded.push_back(make_pair(RunSize{zero_count, bits_count(input[i])}, input[i]));
            zero_count = 0;
        } else {
            zero_count++;
        }
    }
    if (zero_count > 0) {
        encoded.push_back(make_pair(RunSize{0, 0}, 0));
    }
    return encoded;
}


//Hard coded for 8x8
inline vector<int> zigzag_flatten(const vector<vector<int>>& block) {
    vector<int> res(64, 0);
    for (unsigned int y = 0; y < 8; y++)
        for (unsigned int x = 0; x < 8; x++) {
            int index = zigzag.at(y).at(x);
            res[index] = block.at(y).at(x);
        }
    return res;
}

inline void zigzag_inflate(const vector<int> flatten, vector<vector<int>>& block) {
    for (int y = 0; y < 8; y++)
        for (int x = 0; x < 8; x++) {
            int index = zigzag.at(y).at(x);
            block.at(y).at(x) = flatten[index];
        }
}


// B. Bird - 07/02/2020
//Convenience function to wrap around the nasty notation for 2d vectors
template<typename T>
std::vector<std::vector<T> > create_2d_vector(unsigned int outer, unsigned int inner){
    std::vector<std::vector<T> > V {outer, std::vector<T>(inner,T() )};
    return V;
}

//The floating point calculations we use while converting between
//RGB and YCbCr can occasionally yield values slightly out of range
//for an unsigned char (e.g. -1 or 255.9).
//Furthermore, we want to ensure that any conversion uses rounding
//and not truncation (to improve accuracy).
inline unsigned char round_and_clamp_to_char(double v){
    //Round to int
    int i = (int)(v+0.5);
    //Clamp to the range [0,255]
    if (i < 0)
        return 0;
    else if (i > 255)
        return 255;
    return i;
}

/* The exact RGB <-> YCbCr conversion formula used here is the "JPEG style"
   conversion (there is some debate over the best conversion formula)
 */
struct PixelYCbCr;
struct PixelRGB{
    unsigned char r, g, b;
    PixelYCbCr to_ycbcr(); //Implementation is below (since the PixelYCbCr type has to exist before we can fully define this function)
};

struct PixelYCbCr{
    unsigned char Y, Cb, Cr;
    inline PixelRGB to_rgb(){
        return {
            round_and_clamp_to_char(Y + 1.402*(Cr-128.0)),
            round_and_clamp_to_char(Y-0.344136*(Cb-128.0)-0.714136*(Cr-128.0)),
            round_and_clamp_to_char(Y+1.772*(Cb-128.0))
        };
    }
};


inline PixelYCbCr PixelRGB::to_ycbcr(){
    return {
        round_and_clamp_to_char(0.299*r + 0.587*g + 0.114*b),
        round_and_clamp_to_char(128 + -0.168736*r + -0.331264*g + 0.5*b),
        round_and_clamp_to_char(128 + 0.5*r + -0.418688*g + -0.081312*b)
    };
}

#endif
