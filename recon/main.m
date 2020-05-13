
addpath('helper_functions/')  
addpath('models/') 

% Foward Model
psf = im2double(rgb2gray(imread('../RPi_PSF.png')));
psf = imresize(psf, [500, 500], 'box');
psf = psf/norm(psf, 'fro');
[Ny, Nx] = size(psf);

% x = im2double(rgb2gray(imread('small_cal.png')));
x = zeros(Ny, Nx); x(400, 250) = 1;
x = imresize(x, [Ny, Nx], 'box');

y_full = forward_model_2d(x, psf);
% y_full = imnoise(y_full * 1e-12, 'poisson');

mask = im2double(rgb2gray(imread('mask_proto.png')));
mask = imresize(mask, [Ny, Nx], 'nearest');
imshow(mask);
y = y_full.*mask;
%y = im2double(imread('./led_y.png'));
%y = imresize(y, [Ny, Nx], 'nearest');

% subplot(1,3,1); imshow(y_full, []);
subplot(1,3,2); imshow(y, []);
subplot(1,3,3); imshow(mask, []);

% Run inverse model 
opts.fista_iters = 3000;
opts.denoise_method = 'tv'; % options: 'tv', 'non-neg', 'both'
opts.use_gpu = 1;

opts.tv_lambda = 0.005;      % higher: more TV
opts.tv_iters = 100;         % number of inner-loop iterations

% [xout, loss_list] = fista_plain(y, psf, opts);

if opts.use_gpu 
    y = gpuArray(y);
    psf = gpuArray(psf);
    mask = gpuArray(mask);
end


[xout, loss_list] = fista_spectral(y, psf, mask, opts);

subplot(1,3,1); imshow(mask, []);
title('Mask')
subplot(1,3,2); imshow(y, []);
title('Measurement')
subplot(1,3,3); imshow(gather(xout), []);
title('Xout')

imwrite(mat2gray(gather(xout)), 'out.png');
imwrite(mat2gray(gather(y)), 'measurement.png');



function falloff_mask = get_falloff(Ny, Nx, xi, yi)
    falloff_mask = zeros(Ny, Nx);
    for i = 1:4
        f = zeros(Ny, Nx);
        f(round(yi(i)), round(xi(i))) = 1;
        f = imgaussfilt(f, [15, 25], 'FilterSize', [71, 91], 'FilterDomain', 'frequency');
        falloff_mask = falloff_mask + f;
    end
end

