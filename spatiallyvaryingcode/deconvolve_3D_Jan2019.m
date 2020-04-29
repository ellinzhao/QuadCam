%clear all
colors = 'mono'; % not set up to do color yet
ds = 1; % cannot downsample less that 1/2 (dsprior)
ps = 2.9e-3/ds; % pixel size



%%
try
    gpuDevice;
    use_gpu = 1;
catch
    use_gpu = 0;
end

object_close = 0; % change magnification of PSFs (not sure if this works right)
mag = 1; % only does something if object_close=1

estimatebackground = 0; % estimate background fluorescence


disp('Loading principle components...')

savepath = 'example_data/processed_calibration/';
savename = 'example_bilinearweights.mat';
load([savepath savename]);


dsprior = 1; % prior downsampling during registration


pc = pcds;
weights = weightsds;

weights = weights/max(weights(:));
[Ny, Nx, Npc, Nz] = size(pc);

if object_close; Nz = numel(mag); end


disp('warping calibration images...')
if object_close
pcz = zeros(size(pc,1), size(pc,2), Npc, Nz);
    for j = 1:Nz
        for i = 1:Npc
            psf = pc(:,:,i);
            if mag(j) ~= 1
            tform = affine2d([mag(j) 0 0; 0 mag(j) 0; 0 0 1]);
            width = size(psf,2); height = size(psf,1);
            hwarp = imwarp(psf, tform, 'cubic');
            ulx = size(hwarp,2)/2 - width/2;
            uly = size(hwarp,1)/2 - height/2;
            if mag(j) > 1
                psf_warp = imcrop(hwarp, [ulx, uly, width-1, height-1]);
            elseif mag(j) < 1
                pad_size = size(psf)-size(hwarp);
                psf_warp = padarray(hwarp, floor(pad_size/2), 'pre');
                psf_warp = padarray(psf_warp, ceil(pad_size/2), 'post');
            end
            psf = psf_warp;
            end
            pcz(:,:,i,j) = psf;
        end
    end
else
    pcz = pc;
end

%%

disp('downsampling calibration images...')
clear pcds2 weightsds2
for i = 1:Npc
    weightsds2(:,:,i) = imresize(weights(:,:,i),ds*dsprior, 'box');
    for j = 1:Nz
        pcds2(:,:,i,j) = imresize(pcz(:,:,i,j),ds*dsprior, 'box');
    end
end
pcz = pcds2;
weights = weightsds2;
[Ny, Nx, Npc, Nz] = size(pcz);
Ny = Ny/2;
Nx = Nx/2;


pad = @(x,val)padarray(padarray(x, [floor(Ny/2), floor(Nx/2)], val, 'pre'), [ceil(Ny/2), ceil(Nx/2)], val, 'post');
padrep = @(x)padarray(x, [Ny/2, Nx/2], 'replicate', 'both');
cc = floor((Nx/2+1)):floor((3*Nx/2));
rc = floor((Ny/2+1)):floor((3*Ny/2));
crop = @(x)x(rc,cc,:);


% Fourier Transform PCs
disp('Fourier transforming principle components...')
PCz = zeros(size(pcz));

for j = 1:Nz
    pcnorm = norm(pcz(:,:,round(Npc/2),j), 'fro'); % all components are normalized seperately
    for i = 1:Npc
        pp = pcz(:,:,i,j)/pcnorm;
        PCz(:,:,i,j) = fft2(pp);
    end
end

PC_all = PCz;
weights_all = weights;

%% Load raw data

raw_path = 'example_data/raw/';
raw_name = '0426_res03.png';
dctsub_data = 0; % remove low DCT coefficients of raw data

br = double(imread('./b_star_localconv_gauss.png'));

B = fft2(br);

lowfreqmaskb = ones(size(B));
if dctsub_data ~= 0
    lowfreqmaskb(1:dctsub_data/2, 1:dctsub_data/2) = 0;
    lowfreqmaskb(1:dctsub_data/2, end-dctsub_data/2+2:end) = 0;
    lowfreqmaskb(end-dctsub_data/2+2:end, 1:dctsub_data/2) = 0;
    lowfreqmaskb(end-dctsub_data/2+2:end, end-dctsub_data/2+2:end) = 0;
end


br = real(ifft2(B.*lowfreqmaskb));

b = imresize(br, ds, 'box');

%b = b - min(b(:));
%b = max(b,0);
b= b*255/(max(max(b)));
%b= max(0, b-40);
%b = b*255/(2^16-1);
figure(3); imagesc(b); axis image
colorbar

%% Solver inverse problem

if ~estimatebackground

clear PC weights
PC = PC_all;
weights = weights_all;
weights = weights./max(max(sum(weights,3)));
ss = sum(weights,3);
ss(ss < .3) = 1;
for i = 1:size(weights,3)
    weights(:,:,i) = weights(:,:,i)./ss;
end

edgecrop = 10;
Apca = @(x)A_pca_3d(PC, weights, x, pad, crop, use_gpu, edgecrop);
Apcat = @(x)A_adj_pca_3d(PC, weights, x, crop, pad, use_gpu, edgecrop);
if edgecrop~= 0
    bin = b(edgecrop+1:end-edgecrop, edgecrop+1:end-edgecrop);
else
    bin = b;
end
GradErrHandle = @(x)linear_gradient(x, Apca, Apcat, bin);


init_style = 'zero';

switch lower(init_style)
    case('zero')
        xinit = (zeros(Ny*2, Nx*2, Nz));
    case('xhat')
        xinit = xhat;
end

tau = .005;
prox_handle = @(x)soft_nonneg( x, tau );
%prox_handle = @(x)prox_self( x );
h2 = figure(16), clf
h1 = figure(15), clf
options.fighandle = h1;

options.stepsize = 1e-4; %1e-6;
options.convTol = 1e-12;
options.residTol = 50;
options.momentum = 'nesterov';
options.disp_figs = 1;
options.disp_fig_interval = 1;
options.xsize = [Ny Nx];
options.disp_gamma = 1; %1/2.2;
options.disp_crop = @(x)rot90(x(:,:,1),0); %(max(abs(x)/max(x(:)),0)).^options.disp_gamma;
options.known_input = 0;
options.force_real = 1;
options.color_map = 'parula';
options.maxIter = 100;
options.customDisplay = 0;
options.displayFunc = @(x)displayFunc_2plane(rot90(x(:,:,1),2), ...
    rot90(x(:,:,2),2), options.fighandle, ...
    h2, options.color_map);


if use_gpu
    PC = gpuArray(PC);
    weights = gpuArray(weights);
    b = gpuArray(b);
    xinit = gpuArray(xinit);
end
disp('Optimizing...')

[xhat, ~] = proxMin(GradErrHandle, prox_handle, xinit, b, options);

xhat = gather(xhat);

img = (max(abs(xhat)/max(xhat(:)),0));

end


%% Background estimation

if estimatebackground

clear PC weights
PC = PC_all;
weights = weights_all;
weights = weights./max(max(sum(weights,3)));
ss = sum(weights,3);
ss(ss < .3) = 1;
for i = 1:size(weights,3)
    weights(:,:,i) = weights(:,:,i)./ss;
end


edgecrop = 12;
dctest = 10;

Apcabg = @(x, BG) A_pca_3d_bg(PC, weights,x, BG, pad,crop, use_gpu, edgecrop);
Apcabgt = @(x)A_adj_pca_3d_bg(PC, weights, x,dctest, crop, pad,use_gpu, edgecrop);
if edgecrop~= 0
    bin = b(edgecrop+1:end-edgecrop, edgecrop+1:end-edgecrop);
else
    bin = b;
end
%padBG = @(BG) padarray(BG, [2*Ny-dctsub, 2*Nx-dctsub], 'post');
GradErrHandle = @(x)linear_gradient_bg(x, Apcabg, Apcabgt, bin, dctest);


%Apca = @(x)reshape(As*x(:), [Ny, Nx]);
%Apcat = @(x)reshape(As'*x(:), [Ny, Nx]*2);
%GradErrHandle = @(x)linear_gradient(x, Apca, Apcat, b);


init_style = 'zero';

switch lower(init_style)
    case('zero')
        xinit = (zeros(2*Ny, 2*Nx, Nz+1));
    case('xhat')
        xinit = xhat;
end

tau = .05;
prox_handle = @(x)proxNonNeg_BG( x , dctest, Ny, Nx, tau);
h1 = figure(15), clf
h2 = figure(16), clf
options.fighandle = h1;

options.stepsize = 5e-5; %1e-6;
options.convTol = 1e-10;
options.residTol = 10;
options.momentum = 'nesterov';
options.disp_figs = 1;
options.disp_fig_interval = 1;
options.xsize = [Ny Nx];
options.disp_gamma = 1; %1/2.2;
options.disp_crop = @(x)crop(rot90(x(:,:,2),0)); %(max(abs(x)/max(x(:)),0)).^options.disp_gamma;
options.disp_crop = @(x)dctcoeff2bg( x(1:dctsub, 1:dctsub, 2), Ny, Nx ); 
options.known_input = 0;
options.force_real = 1;
options.color_map = 'parula';
options.maxIter = 50;
options.customDisplay = 1;
options.displayFunc = @(x)displayFunc_2plane(rot90(x(:,:,1),0), ...
    dctcoeff2bg( x(1:dctest, 1:dctest, Nz+1), Ny, Nx ), options.fighandle, ...
    h2, options.color_map);


if use_gpu
    PC = gpuArray(PC);
    weights = gpuArray(weights);
    b = gpuArray(b);
    xinit = gpuArray(xinit);
end

disp('Optimizing...')
[xhat, ~] = proxMin(GradErrHandle, prox_handle, xinit, b, options);
xhat = gather(xhat);

% for k = 1:10
%     [xhat, ~] = proxMin(GradErrHandle, prox_handle, xhat, b, options);
% xhat = gather(xhat);
% name = sprintf('bs6_bgest3D_15umpsf0000_11planes_iter%d_ds8.mat', (k+1)*100);
% 
% save(name, 'xhat', 'dctest', 'tau', 'raw_name', 'raw_path', 'mag');
% 
% 
% end

img = max(abs(xhat(:,:,1))/max(max(xhat(:,:,1))),0);

end

        