function [ b ] = A_pca_3d(PC, weights, x, pad, crop, gputrue, edgecrop)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for n = 1:size(PC,4)  % iterate over z-planes
    for m = 1:size(PC,3) % iterate over weight
        if m == 1 && n == 1
            B = PC(:,:,m,n) .* fft2(weights(:,:,m).*x(:,:,n));
        else
            B = B + PC(:,:,m,n) .* fft2(weights(:,:,m).*x(:,:,n));
        end
    end
end
b = crop(fftshift(real(ifft2(B))));
if exist('edgecrop', 'var') || edgecrop ~= 0
    b = b(edgecrop+1:end-edgecrop, edgecrop+1:end-edgecrop);
end


end

