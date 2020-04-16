function Atb = A_adj_pca_3d(PC,weights, x, crop, pad, gputrue, edgecrop)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

padval = mean(x(:));
if exist('edgecrop', 'var') || edgecrop ~= 0
x = padarray(x, [edgecrop edgecrop], padval);
end

if gputrue
    Atb = gpuArray(zeros(size(PC,1), size(PC,2), size(PC,4)));
else
    Atb = zeros(size(PC,1), size(PC,2), size(PC,4));
end

Bp = fft2(pad(x,padval));
for n = 1:size(PC,4)
    for m = 1:size(PC,3)
        PC_conj = conj(PC(:,:,m,n));
        if m == 1
            Atb(:,:,n) = (weights(:,:,m).*fftshift(real(ifft2(PC_conj.*Bp))));
        else
            Atb(:,:,n) = Atb(:,:,n) + (weights(:,:,m).*fftshift(real(ifft2(PC_conj.*Bp))));
        end
    end
end

%Atb = Atb;


end

