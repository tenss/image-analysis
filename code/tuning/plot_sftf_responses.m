stim_params = [Orilist(rand_stim_idx); SFlist(rand_stim_idx); TFlist(rand_stim_idx)]';

stim_frames = 33:16:(32+16*288);

Oris = unique(stim_params(:,1)); 
SFs = unique(stim_params(:,2));
TFs = unique(stim_params(:,3));
nOris = numel(Oris);
nSFs = numel(SFs);
nTFs = numel(TFs);
%%
for indCell = 1:numel(rois)
    resp = reshape(rois(indCell).dfof_corrected(33:end-32), 16, [])';
    [max_resp(indCell), resp_id] = max(mean(resp(:,9:16),2));
    best_stim(indCell,:) = stim_params(resp_id,:);
    stats = regionprops(full(rois(indCell).footprint));
    pos(indCell, :) = stats.Centroid;
end

idx = max_resp>1;

figure; scatter(pos(idx, 1), pos(idx,2), max_resp(idx)*40, ...
    best_stim(idx,3) - best_stim(idx,2), 'filled');
%%
indCell = 99;
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
    imagesc(flipud(squeeze(r(iOri,:,:))'));
    caxis([0 max(r(:))])
end