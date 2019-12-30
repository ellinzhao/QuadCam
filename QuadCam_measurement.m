load 'cam_tforms.mat'
generate_mask = false;

PSF = imread('RPi_PSF.png');
[height, width, channels] = size(PSF);
panorama = zeros([height width 3], 'like', PSF);
blender = vision.AlphaBlender('Operation', 'Binary mask', 'MaskSource', 'Input port');
panoramaView = imref2d([width height]);
mask = true(size(PSF,1), size(PSF,2));
if generate_mask
    PSF = uint8(zeros(size(PSF)));
end
panorama = step(blender, panorama, PSF, mask);

i = 1;
for c = ['A', 'B', 'C', 'D']
    cam = imread(sprintf('2019_11_18-orientation-4000/%s.png', c));
    if generate_mask
        cam = 255*uint8(ones(size(cam)));
    end    
    mask = imwarp(true(size(cam,1),size(cam,2)), cam_tforms(i), 'OutputView', panoramaView);
    cam = imwarp(cam, cam_tforms(i), 'OutputView', panoramaView);
    panorama = step(blender, panorama, cam, mask);
    i = i + 1;
end

imshow(panorama);
filename = 'QuadCam-orientations.png';
if generate_mask
    filename = 'QuadCam-mask.png';
end
imwrite(panorama, filename);