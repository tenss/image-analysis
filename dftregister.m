function xyshifts = dftregister(fft_template, fft_frame, maxshift)
% register a frame using Kuglin and Hines phase correlation method
%
% function xyShifts = dftregister(fft_template, fft_frames, maxShift)
%
%
% Purpose
% Calculate the x/y translation shift between one or more frames and a template frame. 
%
%
% Inputs
% fft_template - The template (fixed) 2-D image to which we will register fft_frames.
% fft_frames - One more more (moving) frames that will be registered to the template
% maxShift - Shifts larger than "maxShift" pixels are not allowed. If empty, no 
%            maximum shift constraint is applied.
% 
% Outputs
% xyShifts - An array of shifts. First row are x shifts and second row y shifts. 
% 
%
% Adapted from Maxime Rio, 2017



% weighting coefficient to balance between phase correlation (alpha = 1)
% and normal correlation (alpha = 0, no normalization)

% adapted from Maxime Rio, 2017

alpha = 0.5;

% compute phase correlation from the normalized cross-spectrum
cs = fft_template .* conj(fft_frame);
cc = ifft2(cs ./ abs(cs).^alpha, 'symmetric');

% constrain maximum shifts found
if ~isempty(maxshift)
    cc(:, maxshift+2:end-maxshift, :) = NaN;
    cc(maxshift+2:end-maxshift, :) = NaN;
end

% split input dimensions
cc_dims = num2cell(size(cc));
[nx, ny] = deal(cc_dims{1:2});

if numel(cc_dims) > 2
    other_dims = cc_dims(3:end);
else
    other_dims = {1};
end

% deduce (x,y)-shift from the maximum in phase correlation
[~, idx] = max(reshape(cc, [], other_dims{:}), [], 1);
[x, y] = ind2sub([nx, ny], idx);

% compensate for circular shifts
x_mask = x > fix(nx / 2);
x(x_mask) = x(x_mask) - nx;

y_mask = y > fix(ny / 2);
y(y_mask) = y(y_mask) - ny;

% compensate for 1-based indexing
x = x - 1;
y = y - 1;

% concatenate (x,y)-shifts
xyshifts = cat(1, x, y);