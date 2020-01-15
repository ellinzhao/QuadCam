cam_tforms = [];
PSF = imread('RPi_PSF.png');

i = 1;
for c = ['A' 'B' 'C' 'D']
    fprintf('Matching image %s', c);
    cam = imread(sprintf('2020_01_14/orientation/%s.png', c));

    [PSF_points, cam_points] = cpselect(cam, PSF, 'Wait', true);
    tform = fitgeotrans(PSF_points, cam_points, 'similarity');
    cam_registered = imwarp(cam, tform, 'OutputView' ,imref2d(size(PSF)));
    
    figure
    imshowpair(cam_registered, PSF, 'blend')
    
    cam_tforms = [cam_tforms tform];
    i = i + 1;
end    

save('cam_tforms_new.mat', 'cam_tforms')