function respmat = aligntrace(trace, events, window)
% align a timeseries to events of interest
%
% trace - timeseries vector
% events - indices of events
% window - vector of timelags around events to extract

idx = (events+min(window))>=1 & (events+max(window))<=length(trace);
events = events(idx);

frameIdx = repmat(window, numel(events), 1) + repmat(events, 1, numel(window));

respmat = trace(frameIdx);

