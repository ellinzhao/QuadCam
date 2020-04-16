
%% 

psf = rgb2gray(imread('./RPi_PSF.png'));
psf = im2double(imresize(psf, [500, 500], 'box'));
psf = psf / norm(psf, 'fro');

num_around_center_y = 3;
num_around_center_x = 5; 

figure, imagesc(psf);
n = 10;
coordinates = zeros(n,2);
hold on
for i=1:n
    [x, y] = ginput(1);
    coordinates(i,:) = [x, y];
    plot(coordinates(:,1), coordinates(:,2), '+');
end
hold off
close 


%% 

[Ny, Nx, Npc, Nz] = size(psfs);
Ny = Ny/2;
Nx = Nx/2;
pad = @(x,val)padarray(padarray(x, [floor(Ny/2), floor(Nx/2)], val, 'pre'), [ceil(Ny/2), ceil(Nx/2)], val, 'post');
cc = floor((Nx/2+1)):floor((3*Nx/2));
rc = floor((Ny/2+1)):floor((3*Ny/2));
crop = @(x)x(rc,cc,:);
edgecrop = 10;
use_gpu = 0;
Apca = @(x)A_pca_3d(psfs, weights_resize, x, pad, crop, use_gpu, edgecrop);


%% Load data

num_around_center_y = 3;
num_around_center_x = 5; 
step = 400;

Nx = 1936; Ny = 1096; % for allocating memory
psfs = zeros(Ny, Nx, (num_around_center_y*2+1),(num_around_center_x*2+1));

for i = 1:(2*(num_around_center_x)+1)
    for j = 1:(2*(num_around_center_y)+1)
        psfs(:,:,j,i) = psf;
    end
end
psf_center = psfs(:,:,num_around_center_y+1, num_around_center_x+1);

%% Downsample

dsprior = 1; % how much to downsample before registration

[Nyorig, Nxorig, ~, ~] = size(psfs);
psfs = imresize(psfs, dsprior, 'box');
psfs = padarray(psfs, [Nyorig/(2/dsprior), Nxorig/(2/dsprior)]);

[Ny, Nx, ~, ~] = size(psfs);
psfc = psfs(:,:,num_around_center_y+1,num_around_center_x+1);


%% Register in concentric circles

mag = 1;
Nz = numel(mag);
Nzcen = 1;

acfunc = @(x)ifftshift(ifft2(abs(fft2(x).^2)));
ccfunc = @(x,y)ifftshift(ifft2(fft2(x).*conj(fft2(y))));

cenindx = num_around_center_x+1;
cenindy = num_around_center_y+1;
psfsr = zeros(size(psfs));
clear xpos ypos


for k = 1:Nz
    psfsr(:,:,cenindy,cenindx, k) = psfs(:,:,num_around_center_y+1,num_around_center_x+1, k);
end
xpos(cenindy, cenindx) = Nx/2;
ypos(cenindy, cenindx) = Ny/2;


for n = 1:max(num_around_center_x, num_around_center_y)
    disp(n)
    num_side_x = min(2*n + 1, 2*num_around_center_x+1);
    num_side_y = min(2*n + 1, 2*num_around_center_y+1);
    
    for i = 1:num_side_x
        for j = 1:num_side_y
            if i == 1 || i == num_side_x || j == 1 || j == num_side_y
                xind = i+max(num_around_center_x-n,0);
                yind = j+max(num_around_center_y-n,0);
                xind = min(max(xind, 1), 2*num_around_center_x+1);
                yind = min(max(yind, 1), 2*num_around_center_y+1);
                
                if i == 1
                    nearestx = xind+1;
                elseif i == num_side_x
                    nearestx = xind-1;
                else
                    nearestx = xind;
                end
                
                if j == 1
                    nearesty = yind+1;
                elseif j == num_side_y
                    nearesty = yind-1;
                else
                    nearesty = yind;
                end
                
                psf_stationary = psfs(:,:,nearesty, nearestx, Nzcen);
                psf_toshift = psfs(:,:,yind,xind, Nzcen);
                cc = ccfunc(psf_stationary, psf_toshift); % first arguement is stationary, second is one to shift
                [m, ind] = max(cc(:));
                [J,I] = ind2sub(size(cc), ind);


                shiftx = I - Nx/2 - (xpos(nearesty, nearestx) - Nx/2);
                shifty = J - Ny/2 - (ypos(nearesty, nearestx) - Ny/2);
                for k = 1:Nz
                temp = circshift(psfs(:,:,yind,xind, k), round([shifty shiftx]*mag(k)/mag(Nzcen)));
                %temp = imtranslate(psfs(:,:,yind,xind, k), round([shiftx shifty]*mag(k)/mag(Nzcen)));
                %temp = circshift(temp, shiftx, 2);
                psfsr(:,:,yind, xind, k) = temp;
                end
                %psfsr(:,:,yind,xind) = imtranslate(psf_toshift, [shiftx, shifty]);
                xpos(yind,xind) = Nx - I + (xpos(nearesty, nearestx) - Nx/2);
                ypos(yind,xind) = Ny - J + (ypos(nearesty, nearestx) - Ny/2);
                distnearx = xpos(nearesty, nearestx) - xpos(yind, xind);
                distneary = ypos(nearesty, nearestx) - ypos(yind, xind);
                
                if xpos(yind,xind) == 0 && ypos(yind,xind)== 0
                    keyboard
                end
                
            end
        end
    end
end


%% Only use subset of indexes after registration

xr = 1:1:2*num_around_center_x+1;
yr = 1:1:2*num_around_center_y+1;

xpos = xpos(yr, xr);
ypos = ypos(yr, xr);
psfsr = psfsr(:,:,yr, xr);
[Ny, Nx, Py, Px, Nz] = size(psfsr);


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


%% Distance from calibration images to diffuser
% "step" should be distance between neighboring samples in um (defined at
% top)

d = 9; % distance diffuser to sensor (mm)
pixsize = 0.0014 / dsprior; % pixels size (mm)

tyy = ypos(3:end,:) - ypos(2:end-1,:);
ty = ypos(3:end,:) - ypos(2:end-1,:);
tyx = xpos(3:end,:) - xpos(2:end-1,:);
tx = xpos(2:end,1:end-1) - xpos(2:end,2:end);
t = mean( [abs(mean(ty(:))), abs(mean(tx(:)))] ); % shift on sensor in pixels
t = t*pixsize;
z = (d/t)*(step*1e-3);
fprintf('\t distance diffuser to source = %3.4f mm\n',  z);


%%

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
weights_all = weights;