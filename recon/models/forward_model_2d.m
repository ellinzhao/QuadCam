function [y] = forward_model_2d(input_image, psf)

[Ny, Nx] = size(psf);   %Get problem size

% Set up pad/crop
p1 = floor(Ny/2);
p2 = floor(Nx/2);
pad2d = @(x)padarray(x,[p1,p2],'both');  %2D padding
crop2d = @(x)x(p1+1:end-p1,p2+1:end-p2,:); %2D cropping

Hs = fftn(ifftshift(pad2d(psf)));  %Compute 3D spectrum
Hfor = @(x)real((ifftn(Hs.*fftn((x)))));

% Run forward model
y = crop2d(Hfor(pad2d(input_image)));


end