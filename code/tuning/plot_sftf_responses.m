load stim_data
% indices of individual grating presentations
% gratings are presented every 16 frames starting from frame 33
stim_frames = 33:16:(32+16*288);
% columns of stim_params are direction, SF, and TF respectively
Oris = unique(stim_params(:,1)); 
SFs = unique(stim_params(:,2));
TFs = unique(stim_params(:,3));
nOris = numel(Oris);
nSFs = numel(SFs);
nTFs = numel(TFs);

%%
% plot an example cell
indCell = 1;
clf;
plot_pos = [1 2 3 6 9 8 7 4];
for iOri = 1:nOris
    for iSF = 1:nSFs
        for iTF = 1:nTFs
            trial_ind = stim_params(:,1)==Oris(iOri) & ...
                        stim_params(:,2)==SFs(iSF) & ...
                        stim_params(:,3)==TFs(iTF);
            r(iOri,iSF,iTF) = mean(rois(indCell).dfof_corrected(...
                stim_frames(trial_ind)+[9:16]));
            rmat(:,iOri,iSF,iTF) = rois(indCell).dfof_corrected(...
                stim_frames(trial_ind)+[1:16]);
        end
    end
    
    subplot(3,3,plot_pos(iOri));
    % plot SF/TF tuning as a colormap, x-axis is SF, y-axis is TF
    imagesc(flipud(squeeze(r(iOri,:,:))'));
    caxis([0 max(r(:))])
end