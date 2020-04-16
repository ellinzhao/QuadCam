%% Load data


%beadcalib = 2; % version 1, and verison 2, (0 for nonbead calibration)


y_center = 2500;
x_center = -25100;


num_around_center_y = 3;
num_around_center_x = 5; 
step = 400;



Nx = 1936; Ny = 1096; % for allocating memory

psfs = zeros(Ny, Nx, (num_around_center_y*2+1),(num_around_center_x*2+1));


% path to data
pathname = 'example_data/calibration/';
filenamebase = 'calib_04-05-19_15umbead_z1200_';

dctsub_psf = 10; % DCT coefficients to remove
camera_baseline = 5; % any baseline above zero from the camera


for i = 1:(2*(num_around_center_x)+1)
    for j = 1:(2*(num_around_center_y)+1)
        xloc_label = x_center+(i-num_around_center_x-1)*step;
        yloc_label = y_center+(j-num_around_center_y-1)*step;

        filename = sprintf([filenamebase 'x_%d_y_%d.png'], xloc_label, yloc_label); 

        
        
        %dx = (i - (num_around_center_x+1))*(step*1e-3)/ps;
        %dy = (j - (num_around_center_y+1))*(step*1e-3)/ps;
        
        input = double(imread([pathname filename]));
        
        if dctsub_psf ~= 0
            dctin = dct2(input);
            dctin(1:dctsub_psf, 1:dctsub_psf) = 0;
            input = idct2(dctin);
            input = max(input,0);
        else
            input = max(input-camera_baseline, 0);
        end
        
             
        psfs(:,:,j,i) = input;
    end
end

psf_center = psfs(:,:,num_around_center_y+1, num_around_center_x+1);
%imagesc(psf_center); axis image

%%

dsprior = 1/2; % how much to downsample before registration

[Nyorig, Nxorig, Py, Px] = size(psfs);
psfs = imresize(psfs, dsprior, 'box');


psfs = padarray(psfs, [Nyorig/(2/dsprior), Nxorig/(2/dsprior)]);


[Ny, Nx, Py, Px] = size(psfs);
psfc = psfs(:,:,num_around_center_y+1,num_around_center_x+1);


%% register in concentric circles

perfectreg = 0;
mag = 1;
Nz = numel(mag);
Nzcen = 1;

%figure(11); clf
acfunc = @(x)ifftshift(ifft2(abs(fft2(x).^2)));
ccfunc = @(x,y)ifftshift(ifft2(fft2(x).*conj(fft2(y))));

cenindx = num_around_center_x+1;
cenindy = num_around_center_y+1;
psfsr = zeros(size(psfs));
clear xpos ypos

if perfectreg
I = Nx/2+1; J = Ny/2+1;
I = round(round(I/delta)*delta);
J = round(round(J/delta)*delta);
shiftx = I - Nx/2;
shifty = J - Ny/2;
psfsr(:,:,cenindy,cenindx) = imtranslate(psfc, [shiftx, shifty]);
xpos(cenindy,cenindx) = Nx - I; ypos(cenindy,cenindx) = Ny - J;
else
    for k = 1:Nz
        psfsr(:,:,cenindy,cenindx, k) = psfs(:,:,num_around_center_y+1,num_around_center_x+1, k);
    end
xpos(cenindy, cenindx) = Nx/2;
ypos(cenindy, cenindx) = Ny/2;
end


for n = 1:max(num_around_center_x, num_around_center_y)
    disp(n)
    num_side_x = min(2*n + 1, 2*num_around_center_x+1);
    num_side_y = min(2*n + 1, 2*num_around_center_y+1);
    %num_ring = num_side*2 + (num_side-2)*2;
    
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

                if perfectreg
                I = round(round((I- (xpos(nearesty, nearestx) - Nx/2))/delta)*delta)+(xpos(nearesty, nearestx) - Nx/2);
                J = round(round((J- (ypos(nearesty, nearestx) - Ny/2))/delta)*delta)+(ypos(nearesty, nearestx) - Ny/2);
                end
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

xpos1= xpos; ypos1 = ypos; psfsr1 = psfsr;

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
%xlim([0 Nx])
%ylim([0 Ny])

%% Distance from calibration images to diffuser
% "step" should be distance between neighboring samples in um (defined at
% top)

d = 3.8; % distance diffuser to sensor (mm)
pixsize = 2.9e-3 / dsprior; % pixels size (mm)

tyy = ypos(3:end,:) - ypos(2:end-1,:);
ty = ypos(3:end,:) - ypos(2:end-1,:);
tyx = xpos(3:end,:) - xpos(2:end-1,:);
tx = xpos(2:end,1:end-1) - xpos(2:end,2:end);
t = mean( [abs(mean(ty(:))), abs(mean(tx(:)))] ); % shift on sensor in pixels
t = t*pixsize;
z = (d/t)*(step*1e-3);
fprintf('\t distance diffuser to source = %3.4f mm\n',  z);


%%
%psfsr(:,:,4,4) = zeros(size(psfsr(:,:,1,1)));

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





%% weights with overlap

setweights_bilinear_interp

disp('Done, saving...')

savepath = 'example_data/processed_calibration/';
savename = 'example_bilinearweights.mat';


save([savepath savename], 'pcds', 'weightsds', 'xpos', 'ypos', 'z', 'd');


disp('saved')
%end
