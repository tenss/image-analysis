function tifpaths = list_tiff_files(stackpath)
    % retrieve paths of TIFF files in a directory

    % list content of the directory
    files = dir(stackpath);

    % find TIFF files
    tif_mask = ~cellfun(@isempty, regexp({files.name}, '\.tiff?$'));
    tifpaths = cellfun(@(t) fullfile(stackpath, t), {files(tif_mask).name}, ...
                       'UniformOutput', false);

    if isempty(tifpaths)
        error('No .tif(f) file found in %s.', stackpath);
    end

    % sort file paths
    [tifpaths, idx_files] = sort(tifpaths);

    % check consistency of file order and dates, to detect problems like bad
    % ordering (e.g. files ending with '_1.tiff' and '_10.tiff')
    [~, idx_date] = sort([files(tif_mask).datenum]);
    if any(idx_files ~= idx_date)
        warning([ ...
            'Files name order and date order are not consistent. ', ...
            'Files might not be loaded in the right order. ', ...
            'See the documentation for more information about this issue.'])
    end
end