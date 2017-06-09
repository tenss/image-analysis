function respmat = aligntrace(trace, events, window)
% align a timeseries to events of interest
%
% trace - timeseries vector
% events - indices of events
% window - vector of timelags around events to extract

frameIdx = repmat(window, numel(events), 1) + repmat(events, 1, numel(window));

respmat = trace(frameIdx);

