function trframe = shiftframe(frame, sx, sy)
% translate a frame by filling an output frame subregion with the input
% frame, which is faster than imstranslate or imdilate-based translation,
% for integer shifts

% adapted from Maxime Rio, 2017

% preallocated result
trframe = zeros(size(frame), 'like', frame);

% fill part of output frame with input frame
[nx, ny, ~] = size(frame);
trframe(max(1, 1+sx):min(nx, nx+sx), max(1, 1+sy):min(ny, ny+sy), :) = ...
    frame(max(1, 1-sx):min(nx-sx, nx), max(1, 1-sy):min(ny-sy, ny), :);
