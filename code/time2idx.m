function idx = time2idx(eventTimestamps, frameTimestamps)

% find the inter-frame interval
ifi = median(diff(frameTimestamps));


nEvents = numel(eventTimestamps);

idx = nan(size(eventTimestamps));

for indE = 1:nEvents
    % do not include events that happened after the end of recording
    if ~isnan(eventTimestamps(indE)) && ...
            eventTimestamps(indE) <= (max(frameTimestamps)+ifi)
        % for each event, find the index of the frame immediately preceeding
        thisEventFrame = find(frameTimestamps<=eventTimestamps(indE), 1, 'last');
        if ~isempty(thisEventFrame)
            idx(indE) = thisEventFrame;
        end
    end
end