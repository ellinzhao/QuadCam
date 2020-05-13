function [X_out,fun_out]=deblur_tv_fista(Bobs,psf,center,lambda,l,u,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA for solving the linear inverse problem with 
% the total variation regularizer and either reflexive or periodic boundary
% conditions
%
% Based on the paper
% Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
% Total Variation Image Denoising and Deblurring Problems"
% -----------------------------------------------------------------------
% Copyright (2008): Amir Beck and Marc Teboulle
% 
% FISTA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% INPUT
%
% Bobs............................. The observed image which is blurred and noisy
% P .................................... PSF of the blurring operator
% center ......................  A vector of length 2 containing the center
%                                           of the PSF
% lambda ...................... Regularization parameter
% l ................................... Lower bound on each of the components
%                                         if no lower bound exists then l=-Inf
% u..................................... Upper bound on each of the components
%                                          if no upper bound exists then u=Inf
% pars.................................Parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.fig ............................... 1 if the image is shown at each
%                                                      iteration, 0 otherwise (Default=1)
% pars.BC .................................. boundary conditions.
%                                                      'reflexive' (default)  or 'periodic
% pars.tv .................................. type of total variation
%                                                      penatly.  'iso' for isotropic (default)
%                                                      and 'l1' for nonisotropic
% pars.mon ............................... 1 if a monotone version of the
%                                                      algorithm is used an 0 otherwise (default)
% pars.denoiseiter .......... number of iterations of the denoising inner
%                                                      problem (default=10)
% OUTPUT
% 
% X_out ......................... Solution of the problem
%                                          min{||A(X)-Bobs||^2+2*lambda*TV(
%                                          X): l <=X_{ij}<=u}
% fun_all .................... Array containing all function values
%                                          obtained in the FISTA method


% Assigning parameres according to pars and/or default values
flag=exist('pars');
if (flag&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if(flag&isfield(pars,'fig'))
    fig=pars.fig;
else
    fig=1;
end
if (flag&isfield(pars,'BC'))
    BC=pars.BC;
else
    BC='reflexive';
end
if (flag&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end
if (flag&isfield(pars,'mon'))
    mon=pars.mon;
else
    mon=0;
end
if (flag&isfield(pars,'denoiseiter'))
    denoiseiter=pars.denoiseiter;
else
    denoiseiter=10;
end

fun_all = [];


[Ny, Nx] = size(psf);   %Get problem size
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

%Bobs = pad2d(Bobs);

switch BC
    case 'diffuser'
        trans = @(X) Hfor(X);
        itrans = @(X) Hadj(X);
        
        Sbig = Hs; 
    otherwise
        error('Invalid boundary conditions should be reflexive or periodic');
end
% computing the two dimensional transform of Bobs

%The Lipschitz constant of the gradient of ||A(X)-Bobs||^2
L=2*max(max(abs(Sbig).^2));


% fixing parameters for the denoising procedure 
clear parsin
parsin.MAXITER=denoiseiter;
parsin.epsilon=1e-5;
parsin.print=0;
parsin.tv=tv;

% initialization
X_iter=Bobs;
Y=X_iter;
t_new=1;

fprintf('***********************************\n');
fprintf('*   Solving with FISTA      **\n');
fprintf('***********************************\n');
fprintf('#iter  fun-val         tv          denoise-iter      relative-dif\n===============================================\n');
for i=1:MAXITER
    % store the old value of the iterate and the t-constant
    X_old=X_iter;
    t_old=t_new;
    % gradient step
    D= pad2d(crop2d(trans(Y)))-Bobs;
    Y=Y-2/L*itrans(D);
    %Y=real(Y);
     
    %invoking the denoising procedure 
    non_neg_only = 'True';
    if non_neg_only == 'True' 
        Z_iter = max(Y, 0);
        iter = 1;
        %tf.nn.relu(vk - self.alpha * grads)
    else
    if (i==1)
        [Z_iter,iter,fun_denoise,P]=denoise_bound_init(Y,2*lambda/L,l,u,[],parsin);
    else
        [Z_iter,iter,fun_denoise]=denoise_bound(Y,2*lambda/L,l,u, parsin); %Xobs,lambda,l,u,pars
    end
    end
    % Compute the total variation and the function value and store it in
    % the function values vector fun_all if exists.
    t=tlv(Z_iter,tv);
    fun_val=norm(trans(Z_iter)-Bobs,'fro')^2+2*lambda*t;

    X_iter=Z_iter;


    fun_all=[fun_all;fun_val];
    
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+t_old/t_new*(Z_iter-X_iter)+(t_old-1)/t_new*(X_iter-X_old);
    
    % printing the information of the current iteration
    fprintf('%3d    %15.5f %15.5f           %3d                  %15.5f\n',i,fun_val,t,iter,norm(X_iter-X_old,'fro')/norm(X_old,'fro'));
    
    if (fig)
        figure(314)
        imshow(X_iter, [])
   
    end
end

fun_out = fun_all;
X_out=X_iter;
