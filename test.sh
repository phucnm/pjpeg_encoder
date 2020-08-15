counter=0
for filename in `find ./test_images -type f`
do
    echo Checking $filename
    ./pjpeg low $filename validate_temp.pjpeg
    ./pjpeg -d validate_temp.pjpeg "${counter}_0low.bmp"
    ./pjpeg medium $filename validate_temp.pjpeg
    ./pjpeg -d validate_temp.pjpeg "${counter}_1medium.bmp"
    ./pjpeg high $filename validate_temp.pjpeg
    ./pjpeg -d validate_temp.pjpeg "${counter}_2high.bmp"
    let counter++
    rm validate_temp.pjpeg
done
