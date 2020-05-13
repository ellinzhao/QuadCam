%% Initialize image properties.
% Change `path` and `save_name` to the appropriate file paths.

path = './test_overlap_small';
save_name = './xcorr_overlap_big.png';
method = 'Canny';

files = imageDatastore(path);
sort_files(files);
n = numel(files.Files);
imtype = readimage(files, 1);
[ny, nx, ~] = size(imtype);
isRGB = size(imtype, 3) == 3;
tforms(n) = affine2d(eye(3));


%% Calculate transformations for image.
tic;
for i = 2:n
    if isRGB
        im1 = rgb2gray(readimage(files, i-1));
        im2 = rgb2gray(readimage(files, i));
    else
        im1 = mat2gray(readimage(files,i-1));
        im2 = mat2gray(readimage(files,i));
    end
    
    if strcmp(method, 'Canny')
        im1 = edge(im1, 'Canny', [], 30);
        im2 = edge(im2, 'Canny', [], 30);
        xcorr = xcorr2_fft(im1, im2);
        [dy, dx] = find(xcorr == max(xcorr(:)), 1);
    else
        im1 = im1 - mean(im1(:));
        im2 = im2 - mean(im2(:));
        ncc = normxcorr2(im2, im1);
        [dy, dx] = find(ncc == max(ncc(:)), 1);
    end
    dy = dy - ny;
    dx = dx - nx;
    tforms(i).T = [1, 0, 0; 0, 1, 0; dx, dy, 1];
end
toc;

for i = 2:numel(tforms)
    tforms(i).T = tforms(i-1).T * tforms(i).T;
end


%% Identify the limits of the image transformations.

xlim = zeros(n, 2);
ylim = zeros(n, 2);
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 nx], [1 ny]);
end

x_min = min([1; xlim(:)]);
x_max = max([nx; xlim(:)]);
y_min = min([1; ylim(:)]);
y_max = max([ny; ylim(:)]);
width  = round(x_max - x_min); 
height = round(y_max - y_min); 

if false % isRGB
    panorama = zeros([height width 3], 'like', imtype);
else
    panorama = zeros([height width], 'like', imtype);
end

blender = vision.AlphaBlender('Operation', 'Binary mask', 'MaskSource', 'Input port');
panoramaView = imref2d([height width], [x_min x_max], [y_min y_max]);


%% Create the panorama.

for i = 1:n
    im = rgb2gray(readimage(files, i));
    warpedImage = imwarp(im, tforms(i), 'OutputView', panoramaView);
    mask = imwarp(true(size(im,1),size(im,2)), tforms(i), 'OutputView', panoramaView);
    panorama = step(blender, panorama, warpedImage, mask);  
end

subplot(1, 3, 1); imshow(im1, []); 
subplot(1, 3, 2); imshow(im2, []); 
subplot(1, 3, 3); imshow(panorama, []); 
imwrite(panorama, save_name);


%% Helper functions.

function im = add_border(im)
    [Ny, Nx] = size(im);
    offset = 100;
    border = 1;
    im(Ny-offset/2:Ny+offset/2) = border;
    im(1:offset) = border;
    im(:, Nx-offset/2:Nx+offset/2) = border;
    im(:, 1:offset) = border;
end


function sort_files(imdatastore)
    file_names = imdatastore.Files;
    str_suffix = regexp(file_names,'\d*','match');
    dub_suffix = str2double(cat(1,str_suffix{:}));
    [~,ii] = sortrows(dub_suffix, size(dub_suffix, 2));
    sorted = file_names(ii);
    imdatastore.Files = sorted;
end

