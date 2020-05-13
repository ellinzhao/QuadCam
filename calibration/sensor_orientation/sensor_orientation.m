%% Preprocessing for transform.

ds = 8;
psf = im2double(rgb2gray(imread('RPi_PSF.png')));
psf = imresize(psf, 1/ds);
psf = psf - 2*mean(psf(:));
psf(psf < 0) = 0;

txs = zeros(4, 1);
tys = zeros(4, 1);
thetas = zeros(4, 1);
i = 1;
for c = ['A' 'B' 'C' 'D']
    fprintf('Matching image %s\n', c);
    sensor = imread(sprintf('led_strip/%s.png', c));
    sensor = im2double(rgb2gray(sensor));
    sensor = imresize(sensor, 1/ds);
    sensor = sensor - 0.5*mean(sensor(:));
    sensor(sensor < 0) = 0;
    [qi, txi, tyi] = get_parameters(psf, sensor);
    thetas(i) = qi;
    txs(i) = txi*ds;
    tys(i) = tyi*ds;
    figure
    warped = imrotate(sensor, qi);
    warped = imtranslate(warped, [txi, tyi], 'OutputView', 'full');
    imshowpair(psf, warped, 'falsecolor');
    i = i + 1;
end   

save('orientation_params.mat', 'txs', 'tys', 'thetas');


%% Helper functions


function [peak, xcorr] = get_ncc_peak(fixed, moving, theta)
    rotated = imrotate(moving, theta);
    xcorr = normxcorr2(rotated, fixed);
    peak = max(xcorr(:));
end

function [q, tx, ty] = get_parameters(fixed, moving)
    tic;
    q = 0;
    [max_corr, max_ncc] = get_ncc_peak(fixed, moving, 0);
    [peak_180, ncc_180] = get_ncc_peak(fixed, moving, 180);
    if (peak_180 > max_corr) 
        max_corr = peak_180;
        max_ncc = ncc_180;
        q = 180;
    end

    offsets = 0:0.25:15;
    for offset = offsets
        % disp(offset);
        [peak_cw, ncc_cw] = get_ncc_peak(fixed, moving, q+offset);
        [peak_ccw, ncc_ccw] = get_ncc_peak(fixed, moving, q-offset);
        if (peak_cw > max_corr) 
            max_corr = peak_cw;
            max_ncc = ncc_cw;
            q = q + offset;
        end
        if (peak_ccw > max_corr)
            max_corr = peak_ccw;
            max_ncc = ncc_ccw;
            q = q - offset;
        end
    end
    toc;

    [ty, tx] = find(max_ncc == max_corr, 1);
    warped = imrotate(moving, q);
    [sy, sx] = size(warped);
    ty = ty - sy; 
    tx = tx - sx;
end

