function [xout, loss_list]=fista_plain(input_image, psf, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements FISTA with non-negativity or with TV 
% Last update: 8/16/2019
%
%
%
% Inputs
% input_image ............. blurred image 
% psf ..................... PSF, same size as blurred image
% opts..................... options file
%  opts.fista_iters ....... Number of FISTA iterations
%  opts.denoise_method .... Either 'non-neg' or 'tv'
%  opts.tv_lambda ......... amount of tv
%  opts.tv_iters .......... number of inner-loop iterations 
%
% Outputs
% xout .................... deblurred image
% loss_list ............... list of losses 


figure(314)

[Ny, Nx] = size(psf);   %Get problem size
L = 100000;
lambda = opts.tv_lambda;

% TV denoising options: 
l = 0;  u = Inf;
clear parsin
parsin.MAXITER=opts.tv_iters;
parsin.epsilon=1e-5;
parsin.print=0;
parsin.tv='iso';


if strcmp(opts.denoise_method, 'tv') == 1
    opts.denoise_method
    prox = @(x)denoise_bound_init(x, 2*lambda/L, l, u, [], parsin);
    loss = @(err, x) norm(err,'fro')^2 + 2*lambda*tlv(x, parsin.tv);
    
elseif strcmp(opts.denoise_method, 'non-neg') == 1
    opts.denoise_method
    prox = @(x)max(x,0);
    loss = @ (err, x) norm(err,'fro')^2;
else
    
end

% Setup convolutional forward op
p1 = floor(Ny/2);
p2 = floor(Nx/2);
pad2d = @(x)padarray(x,[p1,p2],'both');  %2D padding
crop2d = @(x)x(p1+1:end-p1,p2+1:end-p2,:); %2D cropping

vec = @(X)reshape(X,numel(X),1);
Hs = fftn(ifftshift(pad2d(psf)));  %Compute 3D spectrum
Hs_conj = conj(Hs);

Hfor = @(x)real((ifftn(Hs.*fftn((x)))));
Hadj = @(x)real((ifftn(Hs_conj.*fftn((x)))));

padded_input = pad2d(input_image);

%% Start FISTA 
xk = zeros(Ny*2, Nx*2);
vk = zeros(Ny*2, Nx*2);
tk = 1.0;

loss_list = [];

for i=1:opts.fista_iters
    xold = xk;
    vold = vk;
    told = tk;
    
    error = pad2d(crop2d(Hfor(vold))) - padded_input;
    grads = Hadj(error);
    
    xk = prox(vold - 2/L*grads);
    
    tk = 1 + sqrt(1+4*told^2)/2;
    vk = xk + (told-1)/tk *(xk- xold);
    
    loss_i = loss(error, xk);
    %= norm(error,'fro')^2;
    loss_list = [loss_list, loss_i];
    
    subplot(1,2,1)
    imshow(crop2d(xk), [])
    subplot(1,2,2),
    semilogy(loss_list)
    drawnow
    
end

xout = crop2d(xk);

end
