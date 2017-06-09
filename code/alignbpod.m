function trials = alignbpod(datafile)

load(datafile);

frameStarts = [];
lickTimes = cell(SessionData.nTrials,1);
cueTimes = nan(SessionData.nTrials, 2);
rewardTimes = nan(SessionData.nTrials, 1);
punishTimes = nan(SessionData.nTrials, 1);

rewardedTrials = mat2cell([SessionData.TrialTypes == 21]', ...
    ones(SessionData.nTrials,1),1);

for indT = 1:SessionData.nTrials
    thisTrialStart = SessionData.TrialStartTimestamp(indT);
    if isfield(SessionData.RawEvents.Trial{indT}.Events, 'BNC1Low')
        frameStarts = [ frameStarts; ...
            SessionData.RawEvents.Trial{indT}.Events.BNC1Low + ...
                thisTrialStart ];
    end
    
    if isfield(SessionData.RawEvents.Trial{indT}.Events, 'Port4In')
        lickTimes{indT} = SessionData.RawEvents.Trial{indT}.Events.Port4In + ...
                thisTrialStart;
    end
    
    cueTimes(indT,:) = SessionData.RawEvents.Trial{indT}.States.Cue + ...
        thisTrialStart;   
    rewardTimes(indT) = SessionData.RawEvents.Trial{indT}.States.Reward(1) + ...
        thisTrialStart;
    punishTimes(indT) = SessionData.RawEvents.Trial{indT}.States.WhiteNoise(1) + ...
        thisTrialStart;
end

% estimate inter-frame interval
ifi = median(diff(frameStarts));

% while there is at least one skipped frame (with an ifi more than 50%
% above median
while sum(diff(frameStarts) > ifi*1.5)
    % lets find the first skipped frame in the vector
    skippedFramePos = find(diff(frameStarts) > ifi*1.5, 1, 'first') + 1;
    
    % how many frame triggers did we miss?
    skippedFrameNo = ...
        round((frameStarts(skippedFramePos) - frameStarts(skippedFramePos-1)) / ifi) - 1;
    
    % generate timestamps for skipped frames
    skippedFrameTSs = frameStarts(skippedFramePos-1) + [1:skippedFrameNo]' * ifi;
    
    % insert missing timestamps
    frameStarts = [ frameStarts(1:skippedFramePos-1); ...
        skippedFrameTSs; ...
        frameStarts(skippedFramePos:end) ];
end

rewardFrames = mat2cell(time2idx(rewardTimes, frameStarts), ...
    ones(SessionData.nTrials,1),1);
cueFrames = mat2cell(time2idx(cueTimes, frameStarts), ...
    ones(SessionData.nTrials,1),2);
lickFrames = cellfun(@(l) time2idx(l, frameStarts), lickTimes, 'un', 0);
punishFrames = mat2cell(time2idx(punishTimes, frameStarts), ...
    ones(SessionData.nTrials,1),1);

trials = struct('trialType', rewardedTrials, ...
    'cueFrames', cueFrames, ...
    'lickFrames', lickFrames, ...
    'rewardFrames', rewardFrames, ...
    'punishFrames', punishFrames);

