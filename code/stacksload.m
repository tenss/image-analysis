function [stacks, metadata] = stacksload(stackspath, varargin)
    % STACKSLOAD load stacks of images as tensors
    %
    % [stacks, metadata] = stacksload(stackspath, ...)
    %
    % INPUTS
    %   stackspath - stacks paths, as either
    %       1) a file path for a TIFF file
    %       2) a folder path containing images (TIFF files)
    %       3) a file path for binary stack (from IRIS)
    %       4) a cellarray of the previous types
    %
    % NAME-VALUE PAIR INPUTS (optional)
    %   nbplanes - default: []
    %       number of z-planes interleaved in stack frames
    %   nbchannels - default: []
    %       number of channels interleaved in stack frames
    %   forcecell - default: false
    %       return cellarrays, even if one stack path is provided
    %
    % OUTPUTS
    %   stacks - image sequences as [X Y Z Channels Time] array-like objects,
    %       either
    %       1) a TIFFStack object
    %       2) a TensorStack object made of TIFFStack objects
    %       3) a MappedTensor object
    %       4) a cellarray of the previous types (if several stacks)
    %   metadata - (optional) metadata associated with loaded files, as either
    %       1) a [Z Channels Time] structure array of TIFF frames headers
    %       2) an empty array, for IRIS files (no metadata extracted)
    %       3) a cellarray of the previous types (if several stacks)
    %
    % REMARKS
    %   If input folders contain a GRABinfo.mat file from ScanImage, its
    %   content is used if 'nbplanes' and 'nbchannels' values are not provided.
    %
    %   Similarly, ScanImage TIFF file format (v5 and v2016 supported) contains
    %   information that are used if 'nbplanes' and 'nbchannels' values are not
    %   provided.
    %
    %   Returned objects (TIFFStack, TensorStack of TIFFStack or MappedTensor)
    %   are made to access data on demand, from the hard-drive, avoiding
    %   blowing up your computer memory.
    %   They look like normal tensor objects but you *should not* try to
    %   retrieve all frames but get small slices and/or iterate over the
    %   time axis (e.g. to compute a sum of frames).
    %
    % ISSUES
    %   When loading a folder of TIFF files, they are sorted using filenames
    %   before concatenation. This causes problems with filenames that only
    %   differ by trailing numbers.
    %   For example
    %       'file_1.tif' 'file_2.tif' ... 'file_10.tif'
    %   will be sorted as
    %       'file_1.tif' 'file_10.tif' 'file_2.tif' ... 'file_9.tif'
    %   To fix this, you should add zeros in filenames, as follows
    %       'file_01.tif' 'file_02.tif' ... 'file_10.tif'
    %
    % EXAMPLES
    %   % a directory of TIFF files
    %   stack = stacksload('folderpath')
    %
    %   % stack with 4 planes and 2 channels interleaved in frame sequences
    %   stack = stacksload('folderpath', 'nbplanes', 4, 'nbchannels', 2)
    %
    %   % several directories of TIFF files
    %   stacks = stacksload({'folderpath1', 'folderpath2'})
    %
    %   % retrieve stacks and metadata
    %   [stacks, metadata] = stacksload({'folderpath1', 'folderpath2'})
    %
    % SEE ALSO TensorStack, TIFFStack, MappedTensor, stacksmean

    % input checks
    if ~exist('stackspath', 'var')
        error('Missing stackspath argument.');
    end

    if ~iscell(stackspath)
        stackspath = {stackspath};
    end

    if ~iscellstr(stackspath)
        error('Expected stackspath to be a string or a cellarray of strings.');
    elseif any(cellfun(@(x) ~exist(x, 'file'), stackspath))
        error('Some stack paths do not exist.');
    end

    % parse optional inputs
    parser = inputParser;
    posint_attr = {'scalar', 'integer', 'positive'};
    parser.addParameter('nbplanes', [], ...
        @(x) validateattributes(x, {'numeric'}, posint_attr, '', 'nbplanes'));
    parser.addParameter('nbchannels', [], ...
        @(x) validateattributes(x, {'numeric'}, posint_attr, '', 'nbchannels'));
    parser.addParameter('forcecell', false, ...
        @(x) validateattributes(x, {'logical'}, {'scalar'}, '', 'forcecell'));

    parser.parse(varargin{:});
    nbplanes = parser.Results.nbplanes;
    nbchannels = parser.Results.nbchannels;
    forcecell = parser.Results.forcecell;

    % load each stack
    nstacks = numel(stackspath);
    outputs = cell(max(nargout, 1), nstacks);

    for ii=1:nstacks
        % folder of TIFF files
        if isdir(stackspath{ii})
            tifpaths = list_tiff_files(stackspath{ii});
            [outputs{:, ii}] = load_tiff(tifpaths, nbplanes, nbchannels);

        % one TIFF file
        elseif ~isempty(regexp(stackspath{ii}, '\.tiff?', 'once'))
            [outputs{:, ii}] = load_tiff(stackspath(ii), nbplanes, nbchannels);

        % IRIS file
        else
            [outputs{1, ii}] = load_iris(stackspath{ii}, nbplanes, nbchannels);
        end
    end

    % split stack/metadata outputs
    stacks = outputs(1, :);
    if nargout > 1
        metadata = outputs(2, :);
    end

    % return the stack (and its metadata) if there is only one
    if ~forcecell && nstacks == 1
        stacks = stacks{1};
        if nargout > 1
            metadata = metadata{1};
        end
    end
end