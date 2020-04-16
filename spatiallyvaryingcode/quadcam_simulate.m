[Ny, Nx, Npc, Nz] = size(pcds);
Ny = Ny/2;
Nx = Nx/2;
pad = @(x,val)padarray(padarray(x, [floor(Ny/2), floor(Nx/2)], val, 'pre'), [ceil(Ny/2), ceil(Nx/2)], val, 'post');
cc = floor((Nx/2+1)):floor((3*Nx/2));
rc = floor((Ny/2+1)):floor((3*Ny/2));
crop = @(x)x(rc,cc,:);
edgecrop = 0;
use_gpu = 0;
weights = weightsds/max(weightsds(:));

PCz = zeros(size(pcds));
for j = 1:Nz
    pcnorm = norm(pcds(:,:,round(Npc/2),j), 'fro'); % all components are normalized seperately
    for i = 1:Npc
        pp = pcds(:,:,i,j)/pcnorm;
        PCz(:,:,i,j) = fft2(pp);
    end
end
pcds = PCz;

Apca = @(x)A_pca_3d(pcds, weightsds, x, pad, crop, use_gpu, edgecrop);

