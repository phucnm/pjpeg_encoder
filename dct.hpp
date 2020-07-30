//
//  dct.cpp
//  Toy Image Compressor
//
//  Created by Phuc Nguyen on 07/21/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//


#ifndef DCT_HPP
#define DCT_HPP

#include <vector>
#include <array>
#include <math.h>
#include "pjpeg_common.hpp"

using namespace std;

using vvd = vector<vector<double>>;
using vvi = vector<vector<int>>;

// Hardcoded for block 8x8
inline void dct(const vector<vector<unsigned char>>& input, vvd& output) {
    unsigned int i, j, m, n;
    double g;
    double c = 1 / sqrt(2);
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            g = 0;
            double ci = i == 0 ? c : 1;
            double cj = j == 0 ? c : 1;
            for (m = 0; m < 8; m++) {
                for (n = 0; n < 8; n++) {
                    g += input[m][n]
                    * ci * cos(M_PI * (2 * m + 1) * i / 16)
                    * cj * cos(M_PI * (2 * n + 1) * j / 16);
                }
            }
            output[i][j] = g / 4;
        }
    }
}

// Hardcoded for block 8x8
inline void idct(const vvi& input, vector<vector<unsigned char>>& output) {
    unsigned int i, j, m, n;
    double g;
    double c = 1 / sqrt(2);
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            g = 0;
            for (m = 0; m < 8; m++) {
                for (n = 0; n < 8; n++) {
                    double ci = m == 0 ? c : 1;
                    double cj = n == 0 ? c : 1;
                    g += input[m][n]
                    * ci * cos(M_PI * (2 * i + 1) * m / 16)
                    * cj * cos(M_PI * (2 * j + 1) * n / 16);
                }
            }
            output[i][j] = round_and_clamp_to_char(g / 4);
        }
    }
}

inline void quantize(const vvd& input, const vvi& q_vector, vvi& output) {
    for (unsigned int i = 0; i < input.size(); i++) {
        for (unsigned int j = 0; j < input[i].size(); j++) {
            output[i][j] = (int)(input[i][j] / q_vector[i][j]);
        }
    }
}

inline void dequantize(const vvi& input, const vvi& q_vector, vvi& output) {
    for (unsigned int i = 0; i < input.size(); i++) {
        for (unsigned int j = 0; j < input[i].size(); j++) {
            output[i][j] = input[i][j] * q_vector[i][j];
        }
    }
}

#endif
