function stack = load_tiff(tifpaths, nbplanes, nbchannels)
    % load a set of TIFF files, possibly using a GRABinfo.mat file

    % flag to permute z-planes and time axes if it is a Z-stack
    perm_z_t = false;

    % default values for number of z-planes and channels if none found/provided
    if isempty(nbplanes)
        nbplanes = 1;
    end
    if isempty(nbchannels)
        nbchannels = 1;
    end

    % load stacks into one giant tensor
    if numel(tifpaths) > 1
        imgs = cellfun(@(p) TIFFStack(p, false), tifpaths, 'un', false);
        stack = concat_stacks(imgs, nbchannels, nbplanes);

    % or just one TIFFStack
    else
        try
            % try to load TIFF stack, deinterleaving z-planes and channels
            stack = TIFFStack(tifpaths{1}, false, [nbchannels, nbplanes]);
            imgs = {stack};

        catch err
            % re-throw errors non-related to bad deinterleaving
            err_fcn = @(x) strcmp(x.identifier, 'TIFFStack:WrongFrameDims');
            if ~any(cellfun(err_fcn, err.cause))
                rethrow(err);
            end

            % safely reload the TIFF file, compensating the missing frames
            imgs = {TIFFStack(tifpaths{1}, false)};
            stack = concat_stacks(imgs, nbchannels, nbplanes);
        end
    end

    % reorder axes (z-plane and channels and possibly time)
    if ~perm_z_t
        stack = permute(stack, [1, 2, 4, 3, 5]);
    else
        stack = permute(stack, [1, 2, 5, 3, 4]);
    end
end

function tags = getTiffTags(imgs)
    % helper function to retrieve TIFF tags from TIFFStack
    nframes = prod(size(imgs, 3:5));  %#ok<PSIZE>
    tags = getImageTags(imgs, 1:nframes);
end

function stack = concat_stacks(imgs, nbchannels, nbplanes)
    % helper function to safely concatenate and deinterleave TIFFStack objects

    [nx, ny, ~] = size(imgs{1});
    nframes = sum(cellfun(@(x) size(x, 3), imgs));

    % check if deinterleaving will be fine...
    coeff = ceil(nframes / (nbchannels * nbplanes));
    nframes_toadd = coeff * nbchannels * nbplanes - nframes;

    % ... or if an additional dummy stack is needed
    dummy_tensor = {};
    if nframes_toadd > 0
        value = zeros(1, getDataClass(imgs{1}));
        dummy_tensor = ConstantView(value, [nx, ny, nframes_toadd]);
        dummy_tensor = {dummy_tensor};
        warning([ ...
            'Wrong number of frames to deinterleave the stack. ', ...
            'A dummy stack has been added at the end to compensate. ', ...
            'Hence, you should not rely on the content of the last frame.'])
    end

    % concatenate input tensors (and dummy tensor) and re-order dimensions
    stack = TensorStack(3, imgs{:}, dummy_tensor{:});
    stack = reshape(stack, nx, ny, nbchannels, nbplanes, []);
end