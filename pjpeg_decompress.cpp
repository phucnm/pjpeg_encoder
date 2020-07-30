//
//  pjpeg_decompress.cpp
//  Toy Image Compressor
//
//  Created by Phuc Nguyen on 07/21/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include "input_stream.hpp"
#include "bitmap_image.hpp"
#include "pjpeg_common.hpp"
#include "dct.hpp"


void transform_y(const vector<vector<int>>& Y, vector<vector<unsigned char>>& output) {
    int y_h = (int)Y.size();
    int y_w = (int)Y[0].size();
    auto Y_temp = create_2d_vector<int>(y_h, y_w);
    for (unsigned bh = 0; bh < y_h / 8; bh++)
        for (unsigned bw = 0; bw < y_w / 8; bw++) {
            //grab the 8x8 block starting at top left bh, bw
            auto block = create_2d_vector<int>(8, 8);
            for (unsigned y = 0; y < 8; y++)
                for (unsigned x = 0; x < 8; x++) {
                    block.at(y).at(x) = Y.at(bh * 8 + y).at(bw * 8 + x);
                }
            auto q_output = create_2d_vector<int>(8, 8);
            dequantize(block, const_medium_quantize_vector, q_output);
            auto dct_output = create_2d_vector<unsigned char>(8, 8);
            idct(q_output, dct_output);

            for (unsigned h = 0; h < 8; h++)
                for (unsigned w = 0; w < 8; w++) {
                    // crop the paddings
                    int h_idx = bh*8+h;
                    int w_idx = bw*8+w;
                    if (h_idx < output.size() && w_idx < output[0].size()) {
                        Y_temp.at(bh * 8 + h).at(bw * 8 + w) = dct_output.at(h).at(w);
                    }
                }
        }
    for (unsigned y = 0; y < output.size(); y++)
        for (unsigned int x = 0; x < output[0].size(); x++) {
            output.at(y).at(x) = Y_temp.at(y).at(x);
        }
}

//read next 8x8 quantized block
void read_block(InputBitStream &input_stream, vector<vector<int>>& block) {
    //Experiment with huffman decoding
    //First read number of symbols
    int n_syms = (int)input_stream.read_byte();
    unordered_map<int, u32> lengths_map;
    for (unsigned int i = 0; i < n_syms; i++) {
        int symbol = (int)input_stream.read_u32();
        u32 length = (u32)input_stream.read_byte();
        lengths_map[symbol] = length;
    }

    auto code_table = code_lengths_to_code_table(lengths_map);
    // code: symbol
    unordered_map<pair<u16, u16>, int, hash_pair> lut;
    for (auto c: code_table) {
        lut[c.second] = c.first;
    }
    //Start reading bitstream
//    int count = 64;
    int decoded_size = 0;
    for (unsigned int y = 0; y < 8; y++)
        for (unsigned int x = 0; x < 8; x++) {
            //Read 64 coefs
            int code = 0;
            bool found = false;
            int len = 0;
            while (!found) {
                int read = input_stream.read_bit();
                code = code << 1 | read;
                len++;
                auto pair = make_pair(len, code);
                if (lut.find(pair) != lut.end()) {
                    found = true;
                    block.at(y).at(x) = lut[pair];
                    decoded_size += len;
                }
            }
        }
}

// Read symbols, decode Huffman, then perform dequantize and idct
void read_plane(InputBitStream &input_stream, vector<vector<unsigned char>>& plane) {
    //these sizes are padded, to read all input
    int y_h = (int)plane.size();
    int y_w = (int)plane[0].size();
    for (unsigned int bh = 0; bh < y_h / 8; bh++)
        for (unsigned int bw = 0; bw < y_w / 8; bw++) {
            auto block = create_2d_vector<int>(8, 8);
            read_block(input_stream, block);
            auto q_output = create_2d_vector<int>(8, 8);
            dequantize(block, const_medium_quantize_vector, q_output);
            auto dct_output = create_2d_vector<unsigned char>(8, 8);
            idct(q_output, dct_output);

            //put idct block to the plane matrix
            for (unsigned int h = 0; h < 8; h++)
                for (unsigned int w = 0; w < 8; w++) {
                    plane.at(bh * 8 + h).at(bw * 8 + w) = dct_output.at(h).at(w);
                }
        }
}

void transform_and_fill_zzblock(const vector<vector<int>>& quantize_vector,
                                const vector<int>& zz_block,
                                vector<vector<unsigned char>>& plane,
                                int block_h_idx, int block_w_idx) {
    auto block = create_2d_vector<int>(8, 8);
    zigzag_inflate(zz_block, block);
    auto q_output = create_2d_vector<int>(8, 8);
    dequantize(block, quantize_vector, q_output);
    auto dct_output = create_2d_vector<unsigned char>(8, 8);
    idct(q_output, dct_output);

    //put idct block to the plane matrix
    for (unsigned int h = 0; h < 8; h++)
        for (unsigned int w = 0; w < 8; w++) {
            plane.at(block_h_idx * 8 + h).at(block_w_idx * 8 + w) = dct_output.at(h).at(w);
        }
}

void read_plane2(InputBitStream& input_stream, const vector<vector<int>>& quantize_vector, vector<vector<unsigned char>>& plane) {
    int n_syms = (int)input_stream.read_byte();
    unordered_map<RunSize, u32> lengths_map;
    for (unsigned int i = 0; i < n_syms; i++) {
        int run = input_stream.read_byte();
        int size = input_stream.read_byte();
        auto symbol = RunSize{run, size};
        u32 length = (u32)input_stream.read_byte();
        lengths_map[symbol] = length;
    }

    auto code_table = code_lengths_to_code_table(lengths_map);
    // code: symbol
    unordered_map<pair<u16, u16>, RunSize, hash_pair> lut;
    for (auto c: code_table) {
        lut[c.second] = c.first;
    }

    int y_h = (int)plane.size();
    int y_w = (int)plane[0].size();
    for (unsigned int bh = 0; bh < y_h / 8; bh++)
        for (unsigned int bw = 0; bw < y_w / 8; bw++) {
            //after zigzag ordering
            vector<int> zigzag_block(64, 0);
            int num_decoded = 0;
            int count = 0;
            bool eob = false;
            while (!eob) {
                //Read block
                int code = 0;
                bool found = false;
                int len = 0;
                while (!found) {
                    int read = input_stream.read_bit();
                    code = code << 1 | read;
                    len++;
                    auto pair = make_pair(len, code);
                    if (lut.find(pair) != lut.end()) {
                        found = true;
                        RunSize rs = lut[pair];
                        num_decoded += rs.run + 1;
                        if ((rs.run == 0 && rs.size == 0)) {
                            eob = true;
                            transform_and_fill_zzblock(quantize_vector, zigzag_block, plane, bh, bw);
                            break;
                        }
                        //read the diff value
                        int val = input_stream.read_bits(rs.size);
                        int low_positive_bound = 1<<(rs.size - 1);
                        if (val < low_positive_bound) {
                            int min_val_in_bits_length = -((1<<rs.size) - 1);
                            val = min_val_in_bits_length + val;
                        }

                        int run = rs.run;
                        while (--run >= 0) {
                            zigzag_block[count++] = 0;
                        }
                        zigzag_block[count++] = val;
                        if (num_decoded == 64) {
                            eob = true;
                            transform_and_fill_zzblock(quantize_vector, zigzag_block, plane, bh, bw);
                            break;
                        }
                    }
                }
            }

        }
}

void decompress(string input_filename, string output_filename) {
    std::ifstream input_file{input_filename,std::ios::binary};
    InputBitStream input_stream {input_file};

    unsigned int height = input_stream.read_u32();
    unsigned int width = input_stream.read_u32();

    //Read quality code
    vector<vector<int>> quantize_vector = const_medium_quantize_vector;
    int quality_code = input_stream.read_bits(2);
    float quality_factor = 1;
    if (quality_code == LOW_QUALITY_CODE) {
        quality_factor = 2;
    } else if (quality_code == HIGH_QUALITY_CODE) {
        quality_factor = 0.5;
    }
    for (int x = 0; x < 8; x++)
        for (int y = 0; y < 8; y++)
            quantize_vector.at(x).at(y) *= quality_factor;

    int padded_height = next_mul(height, 8);
    int padded_width = next_mul(width, 8);

    int padded_c_h = next_mul((height + 1) / 2, 8);
    int padded_c_w = next_mul((width + 1) / 2, 8);

    auto Y = create_2d_vector<int>(padded_height,padded_width);
    auto Y_plane = create_2d_vector<unsigned char>(padded_height, padded_width);
    auto Cb_plane = create_2d_vector<unsigned char>(padded_c_h, padded_c_w);
    auto Cr_plane = create_2d_vector<unsigned char>(padded_c_h, padded_c_w);
    read_plane2(input_stream, quantize_vector, Y_plane);
    read_plane2(input_stream, quantize_vector, Cb_plane);
    read_plane2(input_stream, quantize_vector, Cr_plane);

    auto imageYCbCr = create_2d_vector<PixelYCbCr>(height,width);
    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            imageYCbCr.at(y).at(x) = {
                Y_plane.at(y).at(x),
                Cb_plane.at(y/2).at(x/2),
                Cr_plane.at(y/2).at(x/2)
            };
        }
    }


    input_stream.flush_to_byte();
    input_file.close();

    bitmap_image output_image {width,height};

    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            auto pixel_rgb = imageYCbCr.at(y).at(x).to_rgb();
            auto [r,g,b] = pixel_rgb;
            output_image.set_pixel(x,y,r,g,b);
        }
    }

    output_image.save_image(output_filename);
}
