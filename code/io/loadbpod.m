function [ trials, frameStarts ] = loadbpod(datafile, nframes)

load(datafile);

lickTimes = cell(SessionData.nTrials,1);
cueTimes = nan(SessionData.nTrials, 2);
rewardTimes = nan(SessionData.nTrials, 1);
punishTimes = nan(SessionData.nTrials, 1);

rewardedTrials = mat2cell([SessionData.TrialTypes == 11]', ...
    ones(SessionData.nTrials,1),1);
frameStarts = [];

for indT = 1:SessionData.nTrials
    thisTrialStart = SessionData.TrialStartTimestamp(indT);
    if isfield(SessionData.RawEvents.Trial{indT}.Events, 'BNC1High')
        frameStarts = [ frameStarts; ...
            SessionData.RawEvents.Trial{indT}.Events.BNC1High + ...
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

trialStarts = mat2cell(SessionData.TrialStartTimestamp', ...
    ones(SessionData.nTrials,1),1);
rewardTimes = mat2cell(rewardTimes, ...
    ones(SessionData.nTrials,1),1);
cueTimes = mat2cell(cueTimes, ...
    ones(SessionData.nTrials,1),2);
punishTimes = mat2cell(punishTimes, ...
    ones(SessionData.nTrials,1),1);

ifis = diff(frameStarts);

frameSkips = find(ifis>median(ifis)*1.5);
for indSkip = 1:numel(frameSkips)
    framesToAdd = ceil(ifis(frameSkips(indSkip)) / median(ifis)) - 1;
    frameStarts = [ frameStarts; ...
        frameStarts(frameSkips(indSkip)) + [1:framesToAdd]'*median(ifis) ];
end

frameStarts = sort(frameStarts);
frameStarts = linspace(min(frameStarts), max(frameStarts), nframes);
trials = struct('trialType', rewardedTrials, ...
    'trialStarts', trialStarts, ...
    'cueTimes', cueTimes, ...
    'lickTimes', lickTimes, ...
    'rewardTimes', rewardTimes, ...
    'punishTimes', punishTimes);

