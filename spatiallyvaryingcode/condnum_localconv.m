% Scene simulation for condition number analysis.

Npoints = [4 16 64 100 225];
delta = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30];
b_scenes = zeros(500, 500, numel(Npoints), numel(delta));
for i = 1:numel(Npoints)
    p = Npoints(i);
    for j = 1:numel(delta)
        scene = get_scene(Ny, Nx, p, delta(j));
        b = Apca(pad(scene, 0));
        b_scenes(:, :, i, j) = b; 
        % results(i, j) = get_cond_num(psf, scene, mask, 0);
    end
end



function A = get_scene(ny, nx, n, d)
    A = zeros(ny, nx);
    ni = floor(sqrt(n));
    nj = n / ni;
    x = floor(.5*(nx - d * (ni - 1)));
    for i = 1:ni
        y = floor(.5*(ny - d * (nj - 1)));
        for j = 1:nj
            A(y, x) = 1.0;
            y = y + d;
        end
        x = x + d;
    end    
end
