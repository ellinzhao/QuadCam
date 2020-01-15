load '2020_01_14/cam_tforms.mat'
mode = 'measurement'; % 'orientation', 'mask', 'measurement'

PSF = imread('RPi_PSF.png');
[height, width, channels] = size(PSF);
panorama = zeros([height width 3], 'like', PSF);
blender = vision.AlphaBlender('Operation', 'Binary mask', 'MaskSource', 'Input port');
panoramaView = imref2d([width height]);
mask = true(size(PSF,1), size(PSF,2));
if ~strcmp(mode, 'orientation')
    PSF = uint8(zeros(size(PSF)));
end
panorama = step(blender, panorama, PSF, mask);

i = 1;
for c = ['A', 'B', 'C', 'D']
    cam = imread(sprintf('2020_01_14/star_shift_more/%s.png', c));
    if strcmp(mode, 'mask')
        cam = 255*uint8(ones(size(cam)));
    end    
    mask = imwarp(true(size(cam,1),size(cam,2)), cam_tforms(i), 'OutputView', panoramaView);
    cam = imwarp(cam, cam_tforms(i), 'OutputView', panoramaView);
    panorama = step(blender, panorama, cam, mask);
    i = i + 1;
end

imshow(panorama);
filename = sprintf('%s.png', mode);
imwrite(panorama, filename);
