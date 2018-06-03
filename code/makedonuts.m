function rois = makedonuts(rois, radius, trim)
% Make donut neuropil ROIs around each cell, excluding labeled cells
% Inputs
%   rois - structure with original ROIs
%   radius - 
cells = zeros(size(rois(1).footprint,1), ...
    size(rois(1).footprint,2));

se = strel('disk', 2);

for indR = 1:size(rois,2)
    cells(:,:) = cells(:,:) + ...
        imdilate(full(rois(1,indR).footprint), se);
end

cells(cells>1) = 1;
if trim>0
    cells(1:trim,:,:) = 1;
    cells(end-trim+1:end,:,:) = 1;
    cells(:,1:trim,:) = 1;
    cells(:,end-trim+1:end,:) = 1;   
end

notCells = 1 - cells;

disk = strel('disk', radius, 8);

for indR = 1:size(rois,2)
    rois(1,indR).donut = imdilate(full(rois(1,indR).footprint), disk) .* ...
        notCells;
    rois(1,indR).donut = sparse(rois(1,indR).donut==1);
end