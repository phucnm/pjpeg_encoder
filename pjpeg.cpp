//
//  main.cpp
//  Toy Image Compressor
//
//  Created by Phuc Nguyen on 07/21/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#include <iostream>
#include "pjpeg_compress.h"
#include "pjpeg_decompress.h"

using namespace std;

int main(int argc, const char * argv[]) {
    if (argc < 4) {
        cerr<<"Usage: "<<argv[0]<<" [-d] [low/medium/high] <input_file> <output_file>"<<endl;
        cerr<<"Example:\nTo compress: "<<argv[0]<<" low input.bmp input.pjpeg"<<endl;
        cerr<<"To decompress: "<<argv[0]<<" input.pjpeg decompressed.bmp"<<endl;
    }

    bool decompress_mode = strcmp(argv[1], "-d") == 0;
    std::string input_filename {argv[2]};
    std::string output_filename {argv[3]};

    if (decompress_mode) {
        decompress(input_filename, output_filename);
    } else {
        std::string quality{argv[1]};
        compress(input_filename, output_filename, quality);
    }

	return 0;
}
