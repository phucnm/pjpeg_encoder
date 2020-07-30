counter=0
for filename in `find ./test_images -type f`
do
    echo Checking $filename
    ./uvg_compress low $filename validate_temp.pjpeg
    ./uvg_decompress validate_temp.uvg "$counter.bmp"
    let counter++
    ./uvg_compress medium $filename  validate_temp.pjpeg
    ./uvg_decompress validate_temp.uvg "$counter.bmp"
    let counter++
    ./uvg_compress high $filename  validate_temp.pjpeg
    ./uvg_decompress validate_temp.uvg "$counter.bmp"
    let counter++
    rm validate_temp.pjpeg
done
