
fixed = rgb2gray(imread('crop.png'));
fixed = imresize(fixed, 0.5);
%fixed = double(edge(fixed, 'Canny', [], 30));
fixed = fixed - mean(fixed(:));
fixed(fixed < 0) = 0;
% fixed = fixed ./ max(fixed(:));

moving = rgb2gray(imread('crop_match.png'));
moving = imresize(moving, 0.5);
moving = moving(1:1000, 1:1000);
%moving = double(edge(moving, 'Canny', [], 30));
moving = moving - mean(moving(:));
moving(moving < 0) = 0;
% moving = moving ./ max(moving(:));

Rfixed = imref2d(size(fixed));
Rfixed.XWorldLimits = [1 size(fixed, 2)];
Rfixed.YWorldLimits = [1 size(fixed, 1)];
Rmoving = imref2d(size(moving));
tformEstimate = imregcorr(moving,fixed,'rigid');
movingReg = imwarp(moving,tformEstimate,'OutputView',Rfixed);
imshowpair(fixed,movingReg,'montage');

optimizer = registration.optimizer.RegularStepGradientDescent;
optimizer.MaximumIterations = 500;
%optimizer.MaximumStepLength = 1e-10;
%optimizer.MinimumStepLength = 1e-20;
%optimizer.GradientMagnitudeTolerance = 10000;
%optimizer.RelaxationFactor = 0.8;
metric = registration.metric.MeanSquares;

[movingRegistered, tform] = imregister(moving, fixed, 'rigid', optimizer, metric, ...
            'PyramidLevels', 1, 'InitialTransformation', tformEstimate, 'DisplayOptimization',true);

%{
ims = imageDatastore({'crop.png', 'crop_match.png'});
psf = rgb2gray(readimage(ims, 1));
points_psf = detectORBFeatures(psf);
[features_psf, points_psf] = extractFeatures(psf, points_psf, 'Method', 'ORB');
n = numel(ims.Files);
tforms(n) = affine2d(eye(3));
im_sizes = zeros(n, 2);

for i = 2:n
    im = rgb2gray(readimage(ims, i));
    im_sizes(i, :) = size(im);
    points = detectORBFeatures(im);
    [features, points] = extractFeatures(im, points, 'Method', 'ORB');
    indexPairs = matchFeatures(features, features_psf, 'Method', 'Exhaustive', 'Unique', true, ...
        'MatchThreshold', 80.0);
    matched = points(indexPairs(:,1), :);
    matched_psf = points_psf(indexPairs(:, 2), :);
    %tforms(i) = estimateGeometricTransform(matched, matched_psf, 'similarity', ...
    %    'Confidence', 99.9, 'MaxNumTrials', 1000);
    tforms(i) = fitgeotrans(matched, matched_psf, 'NonreflectiveSimilarity');
end

% Compute the output limits  for each transform
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 im_sizes(i,2)], [1 im_sizes(i,1)]);
end

im_sizes(1,:) = size(psf);
max_im_size = max(im_sizes);
xMin = min([1; xlim(:)]);
xMax = max([max_im_size(2); xlim(:)]);
yMin = min([1; ylim(:)]);
yMax = max([max_im_size(1); ylim(:)]);

width  = round(xMax - xMin);
height = round(yMax - yMin);
panorama = zeros([height width], 'uint8');

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

% Create the panorama.
for i = 1:n
    im = rgb2gray(readimage(ims, i));
    warpedImage = imwarp(im, tforms(i), 'OutputView', panoramaView);
    mask = imwarp(true(size(im,1),size(im,2)), tforms(i), 'OutputView', panoramaView);
    panorama = step(blender, panorama, warpedImage, mask);
end

% Save the transforms
tforms = tforms(2:n);
save('cam_orientation_Tforms.mat', 'tforms');

% Write and show the camera perspectives/orientations
imwrite(panorama, 'QuadCam_Orienations.png')
figure
imshow(panorama)


function tform = ransac(matched1, matched2, num_trials)
    % RANSAC assuming rigid transformation
    max_inliers = 0;
    best_tform = eye(3);
    
    for i=1:num_trials
        ri = randi(size(matched1, 1));
        rj = randi(size(matched1, 1));
        
    end
end
%}
