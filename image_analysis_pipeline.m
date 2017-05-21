%% first lets load the image stack

% to do this, we will use the TIFFStack class (which you can find on 
% github). TIFFStack provides a convinient way to access multiframe TIFFs
% without having to load the entire file into memory. frames get read from
% disk only when we actually request the data, which is convenient when
% working with stacks that may be many gigabytes in size.

stackpath = 'retinotopy_00002_00001.tif';
stack = TIFFStack(stackpath);

% check the dimensions of our stack
[nx, ny, nt] = size(stack);
fprintf('stack size is: [%d, %d, %d]\n', nx, ny, nt);

%% we are working with small stack here, so we can load it all into memory
stack = stack(:,:,:);

%% lets have a look at a single frame of the stack
figure;
imagesc(stack(:,:,1));
colormap(gray);
axis equal off;

% compare to mean image stack
avgstack = mean(stack, 3);
figure;
imagesc(avgstack);
colormap(gray);
axis equal off;

%% lets look at some raw pixel values
frame = stack(:,:,1);
% plot the distribution of pixel values in a single frame
figure;histogram(frame(:),[-200:10:1000],'normalization','probability');

% it is skewed and there is a sharp peak that corresponds to the image
% offset, which depends on the configuration of the acquisition board and
% PMT amplifiers

% imaging software will often attempt to remove the offset but may or may
% not do it correctly

% let's use a Gaussian mixture model to approximate it, we'll need it later
options = statset('MaxIter',1000);
gmmodel = fitgmdist(double(frame(:)),5,'Options',options);
hold on;plot([-200:1:1000], pdf(gmmodel,[-200:1:1000]')*10, ...
    'LineWidth', 2);

% we'll use the lowest mean gaussian as out estimate of the offset
offset = min(gmmodel.mu);

%% now let's correct for motion artefacts
% cropping the images to avoid the visual stimulation artefacts on the edge
trim = 40;
fft_template = fft2(avgstack(:,trim+1:end-trim));
fft_frame = fft2(double(stack(:,trim+1:end-trim,:)));

xyshifts = dftregister(fft_template, fft_frame, []);

% check the frame shifts
figure;plot(xyshifts');

% correct frame shifts
stack_corrected = zeros(size(stack),'int16');
for ind = 1:nt
    stack_corrected(:,:,ind) = shiftframe(stack(:,:,ind), ...
        xyshifts(1,ind), xyshifts(2,ind));
end
avgreg = mean(stack_corrected,3);

%%
% lets compare the mean images before and after motion correction
figure;
subplot(1,2,1);
imagesc(avgstack), colormap(gray), axis equal off, caxis([offset 1e3]);

subplot(1,2,2);
imagesc(avgreg), colormap(gray), axis equal off, caxis([offset 1e3]);

%%
% select some rois using a simple gui
RoiMaker(avgreg);

%%
% reshape the stack for easy indexing
stack_corrected = reshape(stack_corrected, nx*ny, nt);
% extract activity traces for each roi
for ind = 1:numel(rois)
    mask = reshape(rois(ind).footprint, nx*ny, 1);
    rois(ind).activity = mean(stack_corrected(mask,:)) - offset;
end

%%
% lets look at the activity of a single roi
figure;
subplot(1,2,1); plot(rois(1).activity)

% as with the image offset, we can use a GMM to estimate f0
subplot(1,2,2);
histogram(rois(1).activity,[0:10:500],'normalization','probability');
gmmodel = fitgmdist(rois(1).activity',3,'Options',options);
hold on;plot([0:500], pdf(gmmodel,[0:500]')*10, ...
    'LineWidth', 2);

% loop over rois to calculate dF/F
for ind = 1:numel(rois)
    gmmodel = fitgmdist(rois(ind).activity',3,'Options',options);
    f0 = min(gmmodel.mu);
    rois(ind).dfof = (rois(ind).activity-f0) / f0;
end