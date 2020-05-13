Ny = 500;
Nx = 500;
Sy = 1944;
Sx = 2592;
y = Ny * Sy / 16079;
x = Nx * Sx / 14635;

mask = zeros(Ny, Nx);
for i = 1:4
    f = zeros(Ny, Nx);
    f(Ny/2, Nx/2) = 1;
    f = imgaussfilt(f, [600,700], 'FilterSize', [y - 1, x - 1], 'FilterDomain', 'frequency');
    mask = mask + f;
end
imshow(mask, []);
