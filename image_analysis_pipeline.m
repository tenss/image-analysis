%% first lets load the image stack

% to do this, we will use the TIFFStack class (which you can find on 
% GitHub: https://github.com/DylanMuir/TIFFStack). TIFFStack provides a 
% convenient way to access large multi-frame TIFFs without having to load 
% the entire file into memory. Frames get read from disk only when we 
% actually request the data, which is convenient when working with stacks 
% that may be many gigabytes in size.

stackpath = 'retinotopy_00002_00001.tif';
tsStack = TIFFStack(stackpath);

% check the dimensions of our stack
[nx, ny, nt] = size(tsStack);
fprintf('stack size is: [%d, %d, %d]\n', nx, ny, nt);

% check the class and memory usage of tsStack
whos tsStack
%% we are working with small stack here, so we can load it all into memory
imStack = tsStack(:,:,:);

% check the class and memory usage of imStack
whos imStack

%% lets have a look at a single frame of the stack
figure
imagesc(imStack(:,:,1))
colormap(gray)
axis equal off

% compare to mean image stack
avgStack = mean(imStack, 3);
figure
imagesc(avgStack)
colormap(gray)
axis equal off


%% lets look at some raw pixel values
frame = imStack(:,:,1);
% plot the distribution of pixel values in a single frame
figure; histogram(frame(:), [-200:10:1000], 'normalization', 'probability');

% it is skewed and there is a sharp peak that corresponds to the image
% offset, which depends on the configuration of the acquisition board and
% PMT amplifiers

% imaging software will often attempt to remove the offset but may or may
% not do it correctly

% let's use a Gaussian mixture model to approximate it, we'll need it later
options = statset('MaxIter',1000);
% fit model with 5 Gaussians
gmmodel = fitgmdist(double(frame(:)),5, 'Options', options);
hold on;
% plot the fit to compare with data distribution
plot([-200:1:1000], pdf(gmmodel,[-200:1:1000]')*10, 'LineWidth', 2);

% we'll use the lowest mean Gaussian as our estimate of the offset
offset = min(gmmodel.mu);


%% now let's correct for motion artefacts
% cropping the images to avoid the visual stimulation artefacts on the edge
trimPix = 40;
% the registration algorithm works in the Fourier domain (look up the 
% convolution theorem if you want to know why), so we start by computing
% the 2D Fourier transform of the raw images and registration template
fft_template = fft2(avgStack(:,trimPix+1:end-trimPix));
fft_frame = fft2(double(imStack(:,trimPix+1:end-trimPix,:)));

xyshifts = dftregister(fft_template, fft_frame, []);

% check the frame shifts
figure, plot(xyshifts');

% correct frame shifts
regStack = zeros(size(imStack), 'int16');
for ind = 1:nt
    regStack(:,:,ind) = shiftframe(imStack(:,:,ind), ...
        xyshifts(1,ind), xyshifts(2,ind));
end
regAvg = mean(regStack,3);


%%
% lets compare the mean images before and after motion correction
figure
subplot(1,2,1)
imagesc(avgStack), colormap(gray), axis equal off, caxis([offset 1e3])

subplot(1,2,2)
imagesc(regAvg), colormap(gray), axis equal off, caxis([offset 1e3])

%%
% select some ROIs using a simple GUI
RoiMaker(regAvg);


%%
% reshape the stack for easy indexing
regStack = reshape(regStack, nx*ny, nt);
% extract activity traces for each ROI
for ind = 1:numel(rois)
    mask = reshape(rois(ind).footprint, nx*ny, 1);
    rois(ind).activity = mean(regStack(mask,:)) - offset;
end

% if we accidentally made empty ROIs, get rid of them
emptyIdx = arrayfun(@(r) nnz(r.footprint)==0, rois);
rois = rois(~emptyIdx);
%%
% lets look at the activity of a single ROI
figure
subplot(1,2,1), plot(rois(1).activity)

% as with the image offset, we can use a GMM to estimate f0
subplot(1,2,2)
histogram(rois(1).activity,[0:10:500],'normalization','probability');
% fit a GMM with 3 Gaussians
gmmodel = fitgmdist(rois(1).activity', 3, 'Options', options);
hold on;
plot([0:500], pdf(gmmodel,[0:500]')*10, 'LineWidth', 2);


% loop over ROIs to calculate dF/F
for ind = 1:numel(rois)
    gmmodel = fitgmdist(rois(ind).activity', 3, 'Options', options);
    rois(ind).f0 = min(gmmodel.mu);
    rois(ind).dfof = (rois(ind).activity-rois(ind).f0) / rois(ind).f0;
end

%%
% let's look at the pattern of activity across the population
figure, subplot(2,2,1);
imagesc(cat(1,rois.dfof));
caxis([-1 5])

subplot(2,2,3);
meanAct = mean(cat(1,rois.dfof));
plot(meanAct);

% why does it look stripy?
subplot(2,2,4);
plot(meanAct, rois(3).activity, '.');

% let's attempt to implement a simple neuropil correction 
for ind = 1:numel(rois)
    [~, stats] = robustfit(meanAct, rois(ind).dfof);
    rois(ind).dfof_corrected = stats.resid';
end

% compare to corrected and uncorrected traces
subplot(2,2,2);
imagesc(cat(1,rois.dfof_corrected))
caxis([-1 5])
