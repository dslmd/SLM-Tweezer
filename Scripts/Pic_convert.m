ratio = 125;

before = imread("phase.png"); before = double(before);
after = before * 157.0/255; after = uint8(after); after = im2uint8(after);

imwrite(after,"phase_converted.png");