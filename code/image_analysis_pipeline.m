%% first lets load the image stack

% to do this, we will use the TIFFStack class (which you can find on 
% GitHub: https://github.com/DylanMuir/TIFFStack). TIFFStack provides a 
% convenient way to access large multi-frame TIFFs without having to load 
% the entire file into memory. Frames get read from disk only when we 
% actually request the data, which is convenient when working with stacks 
% that may be many gigabytes in size.

stackpath = 'retinotopy_all.tif';
tsStack = TIFFStack(stackpath);

% check the dimensions of our stack
[nx, ny, nt] = size(tsStack);
fprintf('stack size is: [%d, %d, %d]\n', nx, ny, nt);

% check the class and memory usage of tsStack
whos tsStack

%% if we are working with small stack here, so we can load it all into memory
%imStack = tsStack(:,:,:);
imStack = tsStack;

% check the class and memory usage of imStack
whos imStack

%% lets have a look at a single frame of the stack
figure
subplot(1,2,1)
imagesc(imStack(:,:,1))
colormap(gray)
axis equal off
caxis([0 500] + 32768)
title('single frame')

% compare to mean image stack
avgStack = mean(imStack, 3);
subplot(1,2,2)
imagesc(avgStack)
colormap(gray)
caxis([0 500] + 32768)
axis equal off
title('frame average')

%% lets look at some raw pixel values
frame = imStack(:,:,1);
% plot the distribution of pixel values in a single frame
% bins = [-200:10:1000]; 

% if using retinotopy_all.tif 
bins = [-200:10:1000] + 32768; 

figure, histogram(frame(:), bins, 'normalization', 'probability');

% it is skewed and there is a sharp peak that corresponds to the image
% offset, which depends on the configuration of the acquisition board and
% PMT amplifiers

% imaging software will often attempt to remove the offset but may or may
% not do it correctly

% let's use a Gaussian mixture model to approximate it, we'll need it later
options = statset('MaxIter',1000);
% fit model with 5 Gaussians
gmmodel = fitgmdist(double(frame(:)), 5, 'Options', options);
hold on

% plot the fit to compare with data distribution
% range = [-200:1:1000];
% if using retinotopy_all.tif 
range = [-200:1:1000] + 32768;

plot(range, pdf(gmmodel,range')*10, 'LineWidth', 2);

% we'll use the lowest mean Gaussian as our estimate of the offset
offset = min(gmmodel.mu);
xlabel('gray level')
ylabel('density')

%% now let's correct for motion artefacts
% cropping the images to avoid the visual stimulation artefacts on the edge
trimPix = 40;
% the registration algorithm works in the Fourier domain (look up the 
% convolution theorem if you want to know why), so we start by computing
% the 2D Fourier transform of the raw images and registration template
fft_template = fft2(avgStack(:, trimPix+1:end-trimPix));

xyshifts = zeros(2, nt);
fprintf('starting registration...\n')
for ind = 1:nt
     fft_frame = fft2(double(imStack(:, trimPix+1:end-trimPix, ind)));
     xyshifts(:, ind) = dftregister(fft_template, fft_frame, []);
     if mod(ind,1000)==0
         fprintf('calculating offset of frame %d...\n', ind)
     end
end

% check the frame shifts
figure, plot(xyshifts')
xlabel('frame number')
ylabel('pixel shift')
legend('x shift', 'y shift')

% correct frame shifts
regStack = zeros(size(imStack), 'uint16');
for ind = 1:nt
    regStack(:,:,ind) = shiftframe(imStack(:,:,ind), ...
        xyshifts(1,ind), xyshifts(2,ind));
    if mod(ind,1000)==0
         fprintf('aligning frame %d...\n', ind)
     end
end
regAvg = mean(regStack,3);


%%
% lets compare the mean images before and after motion correction
figure
hAx(1) = subplot(1,2,1);
imagesc(avgStack), colormap(gray), axis equal off, caxis([offset offset+300])

hAx(2) = subplot(1,2,2);
imagesc(regAvg), colormap(gray), axis equal off, caxis([offset offset+300])

linkaxes(hAx);
%%
% select some ROIs using a simple GUI
RoiMaker(regAvg, [offset offset+300]);


%%
% reshape the stack for easy indexing
regStack = reshape(regStack, nx*ny, nt);
% extract activity traces for each ROI
for ind = 1:numel(rois)
    mask = reshape(rois(ind).footprint, nx*ny, 1);
    rois(ind).activity = mean(regStack(mask,:)) - offset;
end

% extra credit:
% if we accidentally made empty ROIs, get rid of them
% emptyIdx = arrayfun(@(r) nnz(r.footprint)==0, rois);
% rois = rois(~emptyIdx);
%%
% lets look at the activity of a single ROI
figure
cellInd = 1;
subplot(1,2,1), plot(rois(cellInd).activity)

% as with the image offset, we can use a GMM to estimate f0
subplot(1,2,2)
histogram(rois(cellInd).activity,[0:10:1000],'normalization','probability');
% fit a GMM with 3 Gaussians
gmmodel = fitgmdist(rois(cellInd).activity', 3, 'Options', options);
hold on;
plot([0:1000], pdf(gmmodel,[0:1000]')*10, 'LineWidth', 2);


% loop over ROIs to calculate dF/F
for ind = 1:numel(rois)
    gmmodel = fitgmdist(rois(ind).activity', 2, 'Options', options);
    rois(ind).f0 = min(gmmodel.mu);
    rois(ind).dfof = (rois(ind).activity-rois(ind).f0) / rois(ind).f0;
end

%%
% let's look at the pattern of activity across the population
figure
imagesc(cat(1,rois.dfof));
caxis([0 5])
xlabel('frames')
ylabel('cells')

%%
% make neuropil ROIs around each cell excluding other labeled cells
rois = makedonuts(rois, 30, 30);
figure
subplot(1,2,1)
imagesc(rois(2).footprint)
title('ROI mask')
axis equal off
colormap gray 

subplot(1,2,2)
imagesc(rois(2).donut)
title('neuropil mask')
axis equal off

% extract neuropil traces for each ROI
for ind = 1:numel(rois)
    mask = reshape(rois(ind).donut, nx*ny, 1);
    rois(ind).neuropil = mean(regStack(mask,:)) - offset;
end

%%
% compare ROI and neuropil fluorescence frame by frame
figure;
subplot(1,2,1);
plot(rois(cellInd).neuropil, rois(cellInd).activity, '.');
[b, stats] = robustfit(rois(cellInd).neuropil, rois(cellInd).activity);
hold on, plot([0 300],[0 300]*b(2) + b(1))
xlim([50 250])
xlabel('neuropil fluorescence')
ylabel('ROI fluorescence')

% let's attempt to implement a simple neuropil correction 
for ind = 1:numel(rois)
    % estimate correction coefficient
    [b, stats] = robustfit(rois(ind).neuropil, rois(ind).activity);
    rois(ind).activity_corrected = stats.resid' + ...
        b(2) * median(rois(ind).neuropil);
    
    % reestimate F0 and dF/F
    gmmodel = fitgmdist(rois(ind).activity_corrected', 2, ...
        'Options', options);
    rois(ind).f0 = min(gmmodel.mu);
    rois(ind).dfof_corrected = (rois(ind).activity_corrected-rois(ind).f0) ...
        / rois(ind).f0;
end

% compare to corrected and uncorrected traces
subplot(1,2,2)
plot(rois(cellInd).dfof)
hold on
plot(rois(cellInd).dfof_corrected)
xlabel('frames')
ylabel('\DeltaF/F')
legend('raw \DeltaF/F', 'corrected \DeltaF/F')

figure
subplot(1,2,1)
imagesc(cat(1,rois.dfof));
caxis([0 5])
xlabel('frames')
ylabel('cells')
title('raw \DeltaF/F')

subplot(1,2,2)
imagesc(cat(1,rois.dfof_corrected))
caxis([0 5])
xlabel('frames')
ylabel('cells')
title('corrected \DeltaF/F')