function trials = loadbpod(datafile)

load(datafile);

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

trialStarts = mat2cell(SessionData.TrialStartTimestamp', ...
    ones(SessionData.nTrials,1),1);
rewardTimes = mat2cell(rewardTimes, ...
    ones(SessionData.nTrials,1),1);
cueTimes = mat2cell(cueTimes, ...
    ones(SessionData.nTrials,1),2);
punishTimes = mat2cell(punishTimes, ...
    ones(SessionData.nTrials,1),1);

trials = struct('trialType', rewardedTrials, ...
    'trialStarts', trialStarts, ...
    'cueTimes', cueTimes, ...
    'lickTimes', lickTimes, ...
    'rewardTimes', rewardTimes, ...
    'punishTimes', punishTimes);
