
%% Load PSF 

psf = rgb2gray(imread('./RPi_PSF.png'));
psf = im2double(imresize(psf, [500, 500], 'box'));
psf = psf / norm(psf, 'fro');

num_around_center_y = 2;
num_around_center_x = 2; 

step = 400; % not used 
Nx = 500; 
Ny = 500;

psfs = zeros(Ny, Nx, (num_around_center_y*2+1), (num_around_center_x*2+1));

for i = 1:(2*(num_around_center_x)+1)
    for j = 1:(2*(num_around_center_y)+1)
        psfs(:, :, j, i) = psf;
    end
end
dsprior = 1; % how much to downsample before registration
[Nyorig, Nxorig, ~, ~] = size(psfs);
psfs = imresize(psfs, dsprior, 'box');
psfs = padarray(psfs, [Nyorig/(2/dsprior), Nxorig/(2/dsprior)]);
psfc = psfs(:,:,num_around_center_y+1,num_around_center_x+1);


%% Select points 
%{
figure, imagesc(psfs(:,:,1));
n = (2*num_around_center_y + 1) * (2*num_around_center_x + 1);
coordinates = zeros(n, 2);
hold on
for i = 1:n
    [x, y] = ginput(1);
    coordinates(i,:) = [x, y];
    plot(coordinates(:,1), coordinates(:,2), '+');
end
hold off
close 

% reshape by row
xpos = reshape(coordinates(:,1).', 2*num_around_center_x + 1, 2*num_around_center_y + 1).';
ypos = reshape(coordinates(:,2).', 2*num_around_center_x + 1, 2*num_around_center_y + 1).';
%}

%% Sample positions

xpos = zeros(2*num_around_center_y + 1, 2*num_around_center_x + 1);
ypos = zeros(2*num_around_center_y + 1, 2*num_around_center_x + 1);

sub_nx = Nx / (2*num_around_center_x + 1);
sub_ny = Ny / (2*num_around_center_y + 1);
for j = 0:(2*num_around_center_y)
    py = Ny/2 + sub_ny/2 + j*sub_ny;
    for i = 0:(2*num_around_center_x)
        px = Nx/2 + sub_nx/2 + i*sub_nx;
        xpos(2*num_around_center_y-j+1, i+1) = px;
        ypos(2*num_around_center_y-j+1, i+1) = py;
    end
end


%% Plot positions of samples

if exist('z_val', 'var'); figure(round(z_val/100)); else
    figure(13); end
 cla; hold on;
for i = 1:size(xpos,1)
    for j = 1:size(xpos,2)
        plot(xpos(i,j), ypos(i,j), 'ro');
    end
end
hold off

axis equal
drawnow


%% Synthetic PSFs

for i = 1:(2*(num_around_center_x)+1)
    for j = 1:(2*(num_around_center_y)+1)
        x = round(xpos(j, i));
        y = round(ypos(j, i));
        ps_shifted = zeros(size(psfs(:, :, j, i)));
        ps_shifted(y, x) = 1;
        psf_shifted = forward_model_2d(ps_shifted, psfs(:, :, j, i));
        % mask = zeros(size(psfs(:, :, j, i)));
        % mask(500, 500) = 1;
        % mask = imgaussfilt(mask, [160, 160], 'FilterSize', [999, 999], 'FilterDomain', 'frequency');
        % mask = mask ./ 6.2394e-06;
        [ny_full, nx_full] = size(psfs(:, :, j, i));
        psf_shifted = psf_shifted; %.* mask;
        psf_shifted(1:250, :) = 0;
        psf_shifted(750:1000, :) = 0;
        psf_shifted(:, 1:250) = 0;
        psf_shifted(:, 750:1000) = 0;
        psfs(:, :, j, i) = psf_shifted;
    end
end


%% Distance from calibration images to diffuser
% "step" should be distance between neighboring samples in um (defined at
% top)

d = 9; % distance diffuser to sensor (mm)
pixsize = 1.4e-3 / dsprior; % pixels size (mm)

tyy = ypos(3:end,:) - ypos(2:end-1,:);
ty = ypos(3:end,:) - ypos(2:end-1,:);
tyx = xpos(3:end,:) - xpos(2:end-1,:);
tx = xpos(2:end,1:end-1) - xpos(2:end,2:end);
t = mean( [abs(mean(ty(:))), abs(mean(tx(:)))] ); % shift on sensor in pixels
t = t*pixsize;
z = (d/t)*(step*1e-3);
fprintf('\t distance diffuser to source = %3.4f mm\n',  z);

%% 
psfsr = psfs;
[Ny, Nx, Py, Px, Nz] = size(psfsr);
psfr_uw = zeros(Ny, Nx, Py*Px, Nz);
weights = zeros(Py, Px, Py*Px);
for i = 1:Py
    %i
    for j = 1:Px
        for k = 1:Nz
        psfr_uw(:,:,j+Px*(i-1), k) = psfsr(:,:,i,j, k);
        end
        %weights(Px-j+1,i,j+Px*(i-1)) = 1;
        weights(i,j,j+Px*(i-1)) = 1;
        
   
    end
end
%pc = psfr_uw;
%sliderDisplayIm(pc)
%sliderDisplayIm(weights)

pc_all = psfr_uw;


%% Weights with overlap

setweights_bilinear_interp

for i = 1:Py*Px
    weightsds(:,:,i) = weightsds(end:-1:1, end:-1:1, i);
end

disp('Done, saving...')
savepath = 'example_data/processed_calibration/';
savename = 'example_bilinearweights.mat';
save([savepath savename], 'pcds', 'weightsds', 'xpos', 'ypos', 'z', 'd');
disp('saved')