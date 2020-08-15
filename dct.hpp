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

static const double S[] = {
    0.353553390593273762200422,
    0.254897789552079584470970,
    0.270598050073098492199862,
    0.300672443467522640271861,
    0.353553390593273762200422,
    0.449988111568207852319255,
    0.653281482438188263928322,
    1.281457723870753089398043,
};

static const double A[] = {
    NAN,
    0.707106781186547524400844,
    0.541196100146196984399723,
    0.707106781186547524400844,
    1.306562964876376527856643,
    0.382683432365089771728460,
};

// DCT type II, scaled. Algorithm by Arai, Agui, Nakajima, 1988.
// See: https://web.stanford.edu/class/ee398a/handouts/lectures/07-TransformCoding.pdf#page=30
//void dct8transform(vector<double> &v) {
//    dct8transform(v.data());
//}

inline void dct8_transform(double vec[]) {
    const double v0 = vec[0] + vec[7];
    const double v1 = vec[1] + vec[6];
    const double v2 = vec[2] + vec[5];
    const double v3 = vec[3] + vec[4];
    const double v4 = vec[3] - vec[4];
    const double v5 = vec[2] - vec[5];
    const double v6 = vec[1] - vec[6];
    const double v7 = vec[0] - vec[7];

    const double v8 = v0 + v3;
    const double v9 = v1 + v2;
    const double v10 = v1 - v2;
    const double v11 = v0 - v3;
    const double v12 = -v4 - v5;
    const double v13 = (v5 + v6) * A[3];
    const double v14 = v6 + v7;

    const double v15 = v8 + v9;
    const double v16 = v8 - v9;
    const double v17 = (v10 + v11) * A[1];
    const double v18 = (v12 + v14) * A[5];

    const double v19 = -v12 * A[2] - v18;
    const double v20 = v14 * A[4] - v18;

    const double v21 = v17 + v11;
    const double v22 = v11 - v17;
    const double v23 = v13 + v7;
    const double v24 = v7 - v13;

    const double v25 = v19 + v24;
    const double v26 = v23 + v20;
    const double v27 = v23 - v20;
    const double v28 = v24 - v19;

    vec[0] = S[0] * v15;
    vec[1] = S[1] * v26;
    vec[2] = S[2] * v21;
    vec[3] = S[3] * v28;
    vec[4] = S[4] * v16;
    vec[5] = S[5] * v25;
    vec[6] = S[6] * v22;
    vec[7] = S[7] * v27;
}

inline void dct8_2d_transform(vector<vector<double>>& v) {
    //row-wise
    for (int i = 0; i < 8; i++) {
        auto data = v.at(i).data();
        dct8_transform(data);
    }
    //column-wise
    for (int j = 0; j < 8; j++) {
        vector<double> c(8);
        for (int i = 0; i < 8; i++) {
            c.at(i) = v.at(i).at(j);
        }
        dct8_transform(c.data());
        for (int i = 0; i < 8; i++) {
            v.at(i).at(j) = c.at(i);
        }
    }
}



//void dct8inverseTransform(vector<double>& v) {
//    dct8inverseTransform(v.data());
//}

// DCT type III, scaled. A straightforward inverse of the forward algorithm.
inline void dct8_inverseTransform(double vec[]) {
    const double v15 = vec[0] / S[0];
    const double v26 = vec[1] / S[1];
    const double v21 = vec[2] / S[2];
    const double v28 = vec[3] / S[3];
    const double v16 = vec[4] / S[4];
    const double v25 = vec[5] / S[5];
    const double v22 = vec[6] / S[6];
    const double v27 = vec[7] / S[7];

    const double v19 = (v25 - v28) / 2;
    const double v20 = (v26 - v27) / 2;
    const double v23 = (v26 + v27) / 2;
    const double v24 = (v25 + v28) / 2;

    const double v7  = (v23 + v24) / 2;
    const double v11 = (v21 + v22) / 2;
    const double v13 = (v23 - v24) / 2;
    const double v17 = (v21 - v22) / 2;

    const double v8 = (v15 + v16) / 2;
    const double v9 = (v15 - v16) / 2;

    const double v18 = (v19 - v20) * A[5];  // Different from original
    const double v12 = (v19 * A[4] - v18) / (A[2] * A[5] - A[2] * A[4] - A[4] * A[5]);
    const double v14 = (v18 - v20 * A[2]) / (A[2] * A[5] - A[2] * A[4] - A[4] * A[5]);

    const double v6 = v14 - v7;
    const double v5 = v13 / A[3] - v6;
    const double v4 = -v5 - v12;
    const double v10 = v17 / A[1] - v11;

    const double v0 = (v8 + v11) / 2;
    const double v1 = (v9 + v10) / 2;
    const double v2 = (v9 - v10) / 2;
    const double v3 = (v8 - v11) / 2;

    vec[0] = (v0 + v7) / 2;
    vec[1] = (v1 + v6) / 2;
    vec[2] = (v2 + v5) / 2;
    vec[3] = (v3 + v4) / 2;
    vec[4] = (v3 - v4) / 2;
    vec[5] = (v2 - v5) / 2;
    vec[6] = (v1 - v6) / 2;
    vec[7] = (v0 - v7) / 2;
}

inline void dct8_2d_inverse_transform(vector<vector<double>>& v, vector<vector<unsigned char>>& o) {
    //row-wise
    for (int i = 0; i < 8; i++) {
        auto data = v.at(i).data();
        dct8_inverseTransform(data);
    }
    //column-wise
    for (int j = 0; j < 8; j++) {
        vector<double> c(8);
        for (int i = 0; i < 8; i++) {
            c.at(i) = v.at(i).at(j);
        }
        dct8_inverseTransform(c.data());
        for (int i = 0; i < 8; i++) {
            v.at(i).at(j) = c.at(i);
        }
    }

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            o.at(i).at(j) = round_and_clamp_to_char(v.at(i).at(j));
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

inline void dequantize(const vvi& input, const vvi& q_vector, vvd& output) {
    for (unsigned int i = 0; i < input.size(); i++) {
        for (unsigned int j = 0; j < input[i].size(); j++) {
            output[i][j] = (double)(input[i][j] * q_vector[i][j]);
        }
    }
}

#endif
