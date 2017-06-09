%% first lets load the image stack

% to do this, we will use the TIFFStack class (which you can find on 
% GitHub: https://github.com/DylanMuir/TIFFStack). TIFFStack provides a 
% convenient way to access large multi-frame TIFFs without having to load 
% the entire file into memory. Frames get read from disk only when we 
% actually request the data, which is convenient when working with stacks 
% that may be many gigabytes in size.

stackpath = 'file_00003.tif';
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
imagesc(imStack(:,:,1))
colormap(gray)
axis equal off

% compare to mean image stack
framesToAverage = 800;
avgStack = mean(imStack(:,:,1:2:framesToAverage*2), 3);
figure
imagesc(avgStack)
colormap(gray)
axis equal off

%% lets look at some raw pixel values
frame = imStack(:,:,1);
% plot the distribution of pixel values in a single frame
bins = [-200:10:500]; 

% if using retinotopy_all.tif 
%bins = [-200:10:1000] + 32768; 

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
range = [-200:1:500];
% if using retinotopy_all.tif 
%range = [-200:1:1000] + 32768;

plot(range, pdf(gmmodel,range')*10, 'LineWidth', 2);

% we'll use the lowest mean Gaussian as our estimate of the offset
offset = 0;%min(gmmodel.mu);


%% now let's correct for motion artefacts
% cropping the images to avoid the visual stimulation artefacts on the edge
trimPix = 40;
% the registration algorithm works in the Fourier domain (look up the 
% convolution theorem if you want to know why), so we start by computing
% the 2D Fourier transform of the raw images and registration template
fft_template = fft2(avgStack(:, trimPix+1:end-trimPix));

xyshifts = zeros(2, nt/2);
for ind = 1:nt/2
     fft_frame = fft2(double(imStack(:, trimPix+1:end-trimPix, ind*2 - 1)));
     xyshifts(:, ind) = dftregister(fft_template, fft_frame, []);
end

% check the frame shifts
figure, plot(xyshifts');

% correct frame shifts
regStack = zeros(size(imStack, 1),size(imStack, 2),framesToAverage,  'int16');
for ind = 1:framesToAverage
    regStack(:,:,ind) = shiftframe(imStack(:,:,ind*2 - 1), ...
        xyshifts(1,ind), xyshifts(2,ind));
end
regAvg = mean(regStack,3);

clear regStack;
%%
% lets compare the mean images before and after motion correction
figure
hAx(1) = subplot(1,2,1);
imagesc(avgStack), colormap(gray), axis equal off, caxis([offset offset+50])

hAx(2) = subplot(1,2,2);
imagesc(regAvg), colormap(gray), axis equal off, caxis([offset offset+50])

linkaxes(hAx);
%%
% select some ROIs using a simple GUI
RoiMaker(regAvg, [offset offset+50]);


%%
% extract fluorescence from rois on chunk of the stack at a time
chunkSize = 400;
nChunks = ceil( nt/(2*chunkSize) );

% vector for storing the second channel, containing bpod sync signals
secondChannel = zeros(nt/2, 1);

for chunk = 1:nChunks
    framesInChunk = ((chunk-1) * chunkSize + 1) : chunk*chunkSize;
    framesInChunk(framesInChunk>nt/2) = [];
    
    % correct movement in a chunk
    regStack = zeros(size(imStack,1), size(imStack,2), numel(framesInChunk), 'uint16');
    for indF = 1:numel(framesInChunk)
        regStack(:,:,indF) = shiftframe(imStack(:,:,framesInChunk(indF)*2 - 1), ...
            xyshifts(1,framesInChunk(indF)), xyshifts(2,framesInChunk(indF)));
        
        secondChannelFrame = imStack(:,:,framesInChunk(indF)*2);
        secondChannel(framesInChunk(indF)) = mean(secondChannelFrame(:));
    end
    
    % reshape the chunk for easy indexing
    regStack = reshape(regStack, nx*ny, numel(framesInChunk));
    
    % extract activity traces for each ROI
    for indR = 1:numel(rois)
        if chunk == 1
            rois(indR).activity = zeros(nt/2, 1);
        end
        
        mask = reshape(rois(indR).footprint, nx*ny, 1);
        rois(indR).activity(framesInChunk) = mean(regStack(mask,:)) - offset;
    end
end

% extra credit:
% if we accidentally made empty ROIs, get rid of them
emptyIdx = arrayfun(@(r) nnz(r.footprint)==0, rois);
rois = rois(~emptyIdx);

%%
% lets look at the activity of a single ROI
figure
cellInd = 1;
subplot(1,2,1), plot(rois(cellInd).activity)

% as with the image offset, we can use a GMM to estimate f0
subplot(1,2,2)
histogram(rois(cellInd).activity,[0:1:200],'normalization','probability');
% fit a GMM with 3 Gaussians
gmmodel = fitgmdist(rois(cellInd).activity, 3, 'Options', options);
hold on;
plot([0:200], pdf(gmmodel,[0:200]')*1, 'LineWidth', 2);

%%
% loop over ROIs to calculate dF/F
for ind = 1:numel(rois)
    gmmodel = fitgmdist(rois(ind).activity, 3, 'Options', options);
    rois(ind).f0 = min(gmmodel.mu);
    rois(ind).dfof = (rois(ind).activity-rois(ind).f0) / rois(ind).f0;
end


%%
% look at trial triggers from the second scanimage channel
figure, plot(secondChannel);

% find threshold crossings corresponding to trial starts
threshold = 1000;
cueStartFrames = find(diff(secondChannel>threshold)==-1) + 1;
% last trial doesn't get recorded in the behavior file
cueStartFrames = cueStartFrames(1:end-1);

% load behavioral data
trials = loadbpod('tenss1_TwoTonePavlovianGoNoGo_Jun09_2017_Session5.mat');

cueTimes = reshape([trials.cueTimes], 2, [])';

% determine the frame rate on the arduino clock
b = regress(cueTimes(:,1), [ones(size(trialStartFrames)) trialStartFrames]);
ifi = b(2);
frameRate = 1/ifi;

% calculate frame numbers of bpod events
for ind = 1:numel(trials)
    trials(ind).cueStartFrame = cueStartFrames(ind);
    trials(ind).lickFrames = cueStartFrames(ind) + ...
        round((trials(ind).lickTimes-trials(ind).cueTimes(1)) / ifi);
    trials(ind).rewardFrame = cueStartFrames(ind) + ...
        round((trials(ind).rewardTimes-trials(ind).cueTimes(1)) / ifi);
    trials(ind).punishFrame = cueStartFrames(ind) + ...
        round((trials(ind).punishTimes-trials(ind).cueTimes(1)) / ifi);
end

%%
% align a cell to licks
lickFrames = cat(1,trials.lickFrames);
cueFrames = [trials.cueStartFrame];

% time window for aligning licks
window = [-19:30];

% only use licks that are at least 10 frames apart - to differentiate lick
% onsets and bouts
lickFrames = lickFrames([true; diff(lickFrames)>10]);

% align dfof response of a the cell to lick events
cellInd = 14;
lickResp = aligntrace(rois(cellInd).dfof, lickFrames, window);

% ta-da!
figure, imagesc(lickResp)