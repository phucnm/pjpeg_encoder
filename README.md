# pjpeg_encoder
This is for experimental purpose, no image viewer.

# Usage 
To compress a bmp (only 24 bit bmp is supported for now due to the bmp lib) to pjpeg,

`./pjpeg <low/medium/high> inputfile.bmp compressed.pjpeg`

To decompress back to a bmp,

`./pjpeg -d compressed.pjpeg decompressed.bmp`

# Demo
The Lena image is compressed in 3 settings: high, medium and low.
![figure](demo/lena.png)
