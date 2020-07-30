//
//  pjpeg_compress.cpp
//  Toy Image Compressor
//
//  Created by Phuc Nguyen on 07/21/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <cstdint>
#include "output_stream.hpp"
#include "bitmap_image.hpp"
#include "pjpeg_common.hpp"
#include "dct.hpp"
#include <unordered_map>
#include <queue>

/*
 Huffman Tree Node
 */
class Node {
public:
    Node(u32 symbol, u32 value) {
        this->symbol = symbol;
        this->value = value;
    }

    u32 symbol;
    u32 value;
    Node* left;
    Node* right;
};

/*
 Huffman Tree builder.
 Provides a canonical codes builder and package-merge builder as static methods.
 */

class Compare
{
public:
     bool operator ()(const Node* lhs, const Node* rhs) const
       {
           return lhs->value > rhs->value;
       }
};

class HuffmanCodingBuilder {
public:
    //This outputs a map with format [symbol: bit length]
    template<typename T>
    static unordered_map<T, u32> package_merge(const unordered_map<T, int>& freqs, int max_bits) {
        vector<pair<vector<T>, u32>> original_freqs;
        for (auto f: freqs) {
            original_freqs.push_back(make_pair(vector<T>{f.first}, f.second));
        }
        if (original_freqs.empty()) {
            return unordered_map<T, u32>{};
        } else if (original_freqs.size() == 1) {
            auto single_sym = original_freqs[0];
            T sym = single_sym.first[0];
            return unordered_map<T, u32>({ {sym, 1} });
        }
//        vector<pair<u32, u32>> res;
        vector<pair<vector<T>, u32>> pkgs(original_freqs);
        sort(pkgs.begin(), pkgs.end(), [](const pair<vector<T>, u32>& a, const pair<vector<T>, u32>& b) -> bool {
            return a.second < b.second;
        });
        u32 num_pick_items = (u32)(2 * original_freqs.size() - 2);
        while (--max_bits) {
            if (pkgs.size() % 2 == 1) {
                pkgs.pop_back();
            }
            vector<pair<vector<T>, u32>> new_pkgs;
            for (size_t i = 0; i < pkgs.size(); i+= 2) {
                vector<T> merged_sym(pkgs[i].first);
                merged_sym.insert(merged_sym.end(), pkgs[i+1].first.begin(), pkgs[i+1].first.end());

                u32 merge_val = pkgs[i].second + pkgs[i+1].second;
                new_pkgs.push_back(make_pair(merged_sym, merge_val));
            }
            pkgs = vector<pair<vector<T>, u32>>(original_freqs);
            pkgs.insert(pkgs.begin(), new_pkgs.begin(), new_pkgs.end());
            sort(pkgs.begin(), pkgs.end(), [](const pair<vector<T>, u32>& a, const pair<vector<T>, u32>& b) -> bool {
                return a.second < b.second;
            });
        }

        unordered_map<T, u32> map;
        for (size_t i = 0; i < num_pick_items; i++) {
            auto v = pkgs[i].first;
            for (T sym: v) {
                if (map.find(sym) == map.end()) {
                    map[sym] = 1;
                } else {
                    map[sym]++;
                }
            }
        }

        return map;
    }
    static vector<u32> generate_from_freqs(vector<u32> freqs) {
        vector<u32> res(freqs.size());

        //build priority queue
        priority_queue<Node*, vector<Node*>, Compare> q;
        for (u32 i = 0; i < freqs.size(); i++) {
            if (freqs[i] > 0) {
                Node* node = new Node(i, freqs[i]);
                node->left = nullptr;
                node->right = nullptr;
                q.push(node);
            }
        }

        Node *root = nullptr;

        if (q.size() == 1) {
            auto node = q.top();
            res[node->symbol] = 1;
            return res;
        }

        while (q.size() > 1) {
            Node *x = q.top();
            q.pop();
            Node * y = q.top();
            q.pop();

            Node* newNode = new Node(-1, x->value + y->value);
            newNode->left = x;
            newNode->right = y;
            root = newNode;
            q.push(newNode);
        }

        tree_traverse(root, res);

        return res;
    }

    static void tree_traverse(Node *root, vector<u32>& result) {
        do_tree_traverse(root, result, 0);
    }

    static void do_tree_traverse(Node *root, vector<u32>& result, u32 code_length) {
        if (root == nullptr) {
            return;
        }
        if (root->left == nullptr && root->right == nullptr) {
            result[root->symbol] = code_length;
            return;
        }
        do_tree_traverse(root->left, result, code_length + 1);
        do_tree_traverse(root->right, result, code_length + 1);
    }
};

//A simple downscaling algorithm using averaging.
std::vector<std::vector<unsigned char> > scale_down(std::vector<std::vector<unsigned char> > source_image, unsigned int source_width, unsigned int source_height, int factor){

    unsigned int scaled_height = (source_height+factor-1)/factor;
    unsigned int scaled_width = (source_width+factor-1)/factor;

    //Note that create_2d_vector automatically initializes the array to all-zero
    auto sums = create_2d_vector<unsigned int>(scaled_height,scaled_width);
    auto counts = create_2d_vector<unsigned int>(scaled_height,scaled_width);

    for(unsigned int y = 0; y < source_height; y++)
        for (unsigned int x = 0; x < source_width; x++){
            sums.at(y/factor).at(x/factor) += source_image.at(y).at(x);
            counts.at(y/factor).at(x/factor)++;
        }

    auto result = create_2d_vector<unsigned char>(scaled_height,scaled_width);
    for(unsigned int y = 0; y < scaled_height; y++)
        for (unsigned int x = 0; x < scaled_width; x++)
            result.at(y).at(x) = (unsigned char)((sums.at(y).at(x)+0.5)/counts.at(y).at(x));
    return result;
}

// Encode huffman a 8x8 quantized block
void output_block(OutputBitStream &output_stream, const vector<vector<int>>& block) {
    unordered_map<int, int> freqs_map;
    for (int y = 0; y < 8; y++)
        for (int x = 0; x < 8; x++) {
            int e = block.at(y).at(x);
            if (freqs_map.find(e) == freqs_map.end()) {
                freqs_map[e] = 1;
            } else {
                freqs_map[e]++;
            }
        }
    // To vector of freqs


    //
    auto lengths_map = HuffmanCodingBuilder().package_merge(freqs_map, 16);
    int est_size = 0;
    int max_len = 0;
    for (auto l: lengths_map) {
        est_size += freqs_map[l.first] * l.second;
        if (l.second > max_len) {
            max_len = l.second;
        }
    }
    est_size /= 8;
    unordered_map<int, pair<u16, u16>> code_table = code_lengths_to_code_table(lengths_map);
    //Since block is 8x8, we are sure num of symbols less than 64 chars, wherese 1 byte = 256 chars.
    output_stream.push_byte((unsigned char)code_table.size());
    //
    for (auto c: code_table) {
        //Symbol
        output_stream.push_bits(c.first, 32);
        // Code length
        output_stream.push_byte(c.second.first);
    }

    int block_encoded_size = 0;
    int z_count = 0;
    for (int y = 0; y < 8; y++)
        for (int x = 0; x < 8; x++) {
            int e = block.at(y).at(x);
            if (e == 0) {
                z_count++;
            }
            auto code_pair = code_table[e];
            output_stream.push_bits(reverse_bits(code_pair.second, code_pair.first), code_pair.first);
            block_encoded_size += code_pair.first;
    }
}

void output_plane(OutputBitStream &output_stream, const vector<vector<int>>& quantize_vector, const vector<vector<unsigned char>>& in) {
    int in_h = (int)in.size();
    int in_w = (int)in[0].size();
    // each block is a vector of pairs of run size and its value
    vector<vector<pair<RunSize, int>>> data;
    for (int bh = 0; bh < in_h / 8; bh++)
        for (int bw = 0; bw < in_w / 8; bw++) {
            //grab the 8x8 block starting at top left bh, bw
            auto block = create_2d_vector<unsigned char>(8, 8);
            for (int y = 0; y < 8; y++)
                for (int x = 0; x < 8; x++) {
                    block.at(y).at(x) = in.at(bh * 8 + y).at(bw * 8 + x);
                }
            auto output = create_2d_vector<double>(8, 8);
            dct(block, output);
            auto q_output = create_2d_vector<int>(8, 8);
            quantize(output, quantize_vector, q_output);

            vector<int> flatten = zigzag_flatten(q_output);
            //Checking
            auto inflate = create_2d_vector<int>(8, 8);
            zigzag_inflate(flatten, inflate);
            assert(inflate == q_output);
            //End checking
            vector<pair<RunSize, int>> rle = rle_encoding(flatten);

            data.push_back(rle);

//            for (int h = 0; h < 8; h++)
//                for (int w = 0; w < 8; w++) {
//                    data.at(bh * 8 + h).at(bw * 8 + w) = q_output.at(h).at(w);
//                }
//            output_block(output_stream, q_output);
        }

    // Count freqs
    unordered_map<RunSize, int> freqs_map;
    for (unsigned int b = 0; b < data.size(); b++) {
        for (auto rs_pair: data[b]) {
            if (freqs_map.find(rs_pair.first) != freqs_map.end()) {
                freqs_map[rs_pair.first]++;
            } else {
                freqs_map[rs_pair.first] = 1;
            }
        }
    }

    auto lengths_map = HuffmanCodingBuilder().package_merge(freqs_map, 16);
    auto code_table = code_lengths_to_code_table(lengths_map);

    output_stream.push_byte((unsigned char)code_table.size());
    //
    for (auto c: code_table) {
        //Symbol
        u8 run = (u8)c.first.run;
        u8 size = (u8)c.first.size;
        output_stream.push_byte(run);
        output_stream.push_byte(size);
        // Code length
        output_stream.push_byte(c.second.first);
    }

    //Encode out data
    for (unsigned int b = 0; b < data.size(); b++) {
        for (auto rs_pair: data[b]) {
            auto code_pair = code_table[rs_pair.first];
            //output code
            output_stream.push_bits(reverse_bits(code_pair.second, code_pair.first), code_pair.first);
            //output value
            //e.g. min of 3 bits is -7
            // min of 5 bits is -31
            int min_val_in_bits_length = -((1<<rs_pair.first.size) - 1);
            if (rs_pair.second > 0) {
                output_stream.push_bits(rs_pair.second, rs_pair.first.size);
            } else {
                // if num is negative, we output the diff
                int diff = rs_pair.second - min_val_in_bits_length;
                output_stream.push_bits(diff, rs_pair.first.size);
            }
        }
    }
}

void output(OutputBitStream& output_stream,
            const vector<vector<int>>& quantize_vector,
            const vector<vector<unsigned char>>&Y,
            const vector<vector<unsigned char>>&Cb,
            const vector<vector<unsigned char>>&Cr) {
    output_plane(output_stream, quantize_vector, Y);
    output_plane(output_stream, quantize_vector, Cb);
    output_plane(output_stream, quantize_vector, Cr);
}

void compress(string input_filename, string output_filename, string quality) {
    bitmap_image input_image {input_filename};

    unsigned int height = input_image.height();
    unsigned int width = input_image.width();

    //Read the entire image into a 2d array of PixelRGB objects
    //(Notice that height is the outer dimension, so the pixel at coordinates (x,y)
    // must be accessed as imageRGB.at(y).at(x)).
    std::vector<std::vector<PixelYCbCr>> imageYCbCr = create_2d_vector<PixelYCbCr>(height,width);


    for(unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            auto [r,g,b] = input_image.get_pixel(x,y);
            PixelRGB rgb_pixel {r,g,b};
            imageYCbCr.at(y).at(x) = rgb_pixel.to_ycbcr();
        }
    }

    std::ofstream output_file{output_filename,std::ios::binary};
    OutputBitStream output_stream {output_file};

    //Placeholder: Use a simple bitstream containing the height/width (in 32 bits each)
    //followed by the entire set of values in each colour plane (in row major order).

    output_stream.push_u32(height);
    output_stream.push_u32(width);

    //Extract the Y values
    int padded_height = next_mul(height, 8);
    int padded_width = next_mul(width, 8);
    auto Y = create_2d_vector<unsigned char>(padded_height, padded_width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Y.at(y).at(x) = imageYCbCr.at(y).at(x).Y;
    //            output_stream.push_byte(imageYCbCr.at(y).at(x).Y);

    //Extract the Cb plane into its own array
    auto Cb = create_2d_vector<unsigned char>(height,width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Cb.at(y).at(x) = imageYCbCr.at(y).at(x).Cb;
    auto Cb_scaled = scale_down(Cb,width,height,2);
    int cb_padded_h = next_mul((int)Cb_scaled.size(), 8);
    int cb_padded_w = next_mul((int)Cb_scaled[0].size(), 8);
    auto Cb_padded = create_2d_vector<unsigned char>(cb_padded_h, cb_padded_w);
    for(unsigned int y = 0; y < Cb_scaled.size(); y++)
        for (unsigned int x = 0; x < Cb_scaled[0].size(); x++)
            Cb_padded.at(y).at(x) = Cb_scaled.at(y).at(x);

    //Extract the Cr plane into its own array
    auto Cr = create_2d_vector<unsigned char>(height,width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Cr.at(y).at(x) = imageYCbCr.at(y).at(x).Cr;
    auto Cr_scaled = scale_down(Cr,width,height,2);
    int cr_padded_h = next_mul((int)Cr_scaled.size(), 8);
    int cr_padded_w = next_mul((int)Cr_scaled[0].size(), 8);
    auto Cr_padded = create_2d_vector<unsigned char>(cr_padded_h, cr_padded_w);
    for(unsigned int y = 0; y < Cr_scaled.size(); y++)
        for (unsigned int x = 0; x < Cr_scaled[0].size(); x++)
            Cr_padded.at(y).at(x) = Cr_scaled.at(y).at(x);

    //Prepare quantization vector
    // For low just /2, high * 2 otherwise leave it as is.
    vector<vector<int>> quantize_vector = const_medium_quantize_vector;
    int quality_code = MEDIUM_QUALITY_CODE;
    float quality_factor = 1;
    if (quality.compare("low") == 0) {
        quality_factor = 2;
        quality_code = LOW_QUALITY_CODE;
    } else if (quality.compare("high") == 0) {
        quality_factor = 0.5;
        quality_code = HIGH_QUALITY_CODE;
    }
    for (unsigned int x = 0; x < 8; x++)
        for (unsigned int y = 0; y < 8; y++)
            quantize_vector.at(x).at(y) *= quality_factor;

    //Output the quality into bitstream
    output_stream.push_bits(quality_code, 2);

    //Output the 3 planes
    output(output_stream, quantize_vector, Y, Cb_padded, Cr_padded);

    output_stream.flush_to_byte();
    output_file.close();
}
