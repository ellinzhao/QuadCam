load 'orientation_params.mat'

psf = rgb2gray(imread('RPi_PSF.png'));
[height, width] = size(psf);
panorama = double(zeros([height width]));
blender = vision.AlphaBlender('Operation', 'Binary mask', 'MaskSource', 'Input port');
panoramaView = imref2d([width height]);
mask = true(size(psf,1), size(psf,2));
psf = zeros(size(psf));
panorama = step(blender, panorama, psf, mask);

i = 1;
for c = ['A', 'B', 'C', 'D']
    sensor = imread(sprintf('led_plus/%sm.png', c));
    sensor = im2double(rgb2gray(sensor));
    q = thetas(i);
    tx = txs(i);
    ty = tys(i);
    mask = warp(true(size(sensor)), tx, ty, q);
    panorama = step(blender, panorama, double(mask), mask);
    i = i + 1;
end

imshow(panorama);
imwrite(panorama, 'mask.png');


function warped = warp(im, tx, ty, q)
    warped = imrotate(im, q);
    warped = imtranslate(warped, [tx, ty], 'OutputView', 'full');
end
