%%
% load behavioral data
[ trials, frameStarts] = ...
    loadbpod('behavior/Dummy Subject_TwoOdors_May09_2018_Session15.mat', nt);

cueTimes = reshape([trials.cueTimes], 2, [])';

% determine the frame rate on the arduino clock
ifi = median(diff(frameStarts));
frameRate = 1/ifi;

% calculate frame numbers of bpod events
for ind = 1:numel(trials)
    trials(ind).cueStartFrame = find(...
        trials(ind).cueTimes(1)-frameStarts<ifi & trials(ind).cueTimes(1)-frameStarts>=0,...
        1,'first');
    if ~isempty(trials(ind).lickTimes)
        trials(ind).lickFrames = ...
            cell2mat(arrayfun(@(t) find(...
            t-frameStarts<ifi & t-frameStarts>=0,1,'first'), ...
            trials(ind).lickTimes,'un',0));
    end
    if ~isempty(trials(ind).rewardTimes)
        trials(ind).rewardFrame = ...
            find(...
            trials(ind).rewardTimes(1)-frameStarts<ifi & ...
            trials(ind).rewardTimes(1)-frameStarts>=0,1,'first');
    end
    if ~isempty(trials(ind).punishTimes) && ~isnan(trials(ind).punishTimes)
        trials(ind).punishFrame = ...
            find(trials(ind).punishTimes(1)-frameStarts<ifi & ...
            trials(ind).punishTimes(1)-frameStarts>=0,1,'first');
    end
end

%%
% align a cell to licks
%lickFrames = cat(1,trials.lickFrames);
rewardFrames = [trials.rewardFrame]';
rewardFrames = rewardFrames(~isnan(rewardFrames));
cueStartFrames = [trials.cueStartFrame];
unrewardedCues = [trials([trials.trialType]==0).cueStartFrame];
rewardedCues = [trials([trials.trialType]==1).cueStartFrame];

lickTrials = arrayfun(@(t) ~isempty(t.lickTimes) && min(t.lickTimes-t.cueTimes(1))<3, trials);

rewardedCuesLick = [trials([trials.trialType]==1 & lickTrials').cueStartFrame];
rewardedCuesNolick = [trials([trials.trialType]==1 & ~lickTrials').cueStartFrame];
unrewardedCuesLick = [trials([trials.trialType]==0 & lickTrials').cueStartFrame];
unrewardedCuesNolick = [trials([trials.trialType]==0 & ~lickTrials').cueStartFrame];
% time window for aligning licks
window = [-19:40];

% only use licks that are at least 10 frames apart - to differentiate lick
% onsets and bouts
lickFrames = lickFrames([true; diff(lickFrames)>10]);

figure

%%
% align dfof response of a the cell to lick events
for cellInd = 1:numel(rois)
    
clim = [0 3];
lickResp = aligntrace(rois(cellInd).dfof_corrected, lickFrames, window);

rewardResp = aligntrace(rois(cellInd).dfof_corrected, rewardFrames, window);


rewardedCueResp = aligntrace(rois(cellInd).dfof_corrected, rewardedCues', window);
rewardedCueLickResp = aligntrace(rois(cellInd).dfof_corrected, rewardedCuesLick', window);
rewardedCueNolickResp = aligntrace(rois(cellInd).dfof_corrected, rewardedCuesNolick', window);

unrewardedCueResp = aligntrace(rois(cellInd).dfof_corrected, unrewardedCues', window);
unrewardedCueLickResp = aligntrace(rois(cellInd).dfof_corrected, unrewardedCuesLick', window);
unrewardedCueNolickResp = aligntrace(rois(cellInd).dfof_corrected, unrewardedCuesNolick', window);

% ta-da!
subplot(4,1,1), imagesc(unrewardedCueResp)
%hold on, line(repmat(find(window==0),2,1), [0 numel(unrewardedCues)], 'Color', 'w');
set(gca,'XTick',[20:20:60], 'XTickLabel', num2str(window(20:20:60)'*ifi,2), ...
    'TickDir', 'out');
xlabel('Time from odor A')
caxis(clim)

subplot(4,1,2), imagesc(rewardedCueResp)
%hold on, line(repmat(find(window==0),2,1), [0 numel(unrewardedCues)], 'Color', 'w');
set(gca,'XTick',[20:20:60], 'XTickLabel', num2str(window(20:20:60)'*ifi,2), ...
    'TickDir', 'out');
xlabel('Time from tone A')
caxis(clim)

subplot(4,1,3), imagesc(lickResp)
%hold on, line(repmat(find(window==0),2,1), [0 numel(rewardedCues)], 'Color', 'w');
set(gca,'XTick',[20:20:60], 'XTickLabel', num2str(window(20:20:60)'*ifi,2), ...
    'TickDir', 'out');
xlabel('Time from tone B');
caxis(clim)

keyboard;
end