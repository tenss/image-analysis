function trframe = shiftframe(frame, sx, sy)
% Fast frame translation
% 
% function trFrame = shiftframe(frame, sx, sy)
%
% 
% Purpose
% Translate a frame by filling an output frame subregion with the input
% frame, which is faster than imstranslate or imdilate-based translation,
% for integer shifts
%
%
% Inputs
% frame - the frame to shift
% sx - x shift
% sy - y shift
%
%
% Outputs
% trFrame - the translated frame. 
%
%
% Adapted from Maxime Rio, 2017

% preallocated result
trframe = zeros(size(frame), 'like', frame);

% fill part of output frame with input frame
[nx, ny, ~] = size(frame);
trframe(max(1, 1+sx):min(nx, nx+sx), max(1, 1+sy):min(ny, ny+sy), :) = ...
    frame(max(1, 1-sx):min(nx-sx, nx), max(1, 1-sy):min(ny-sy, ny), :);
