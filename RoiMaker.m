classdef RoiMaker < handle


    % A simple class that makes a GUI to interactively place elliptical ROIs on an image
    %
    %
    % Run using: 
    % >>  RoiMaker(2D_imageData);
    %
    % e.g.
    % >>  RoiMaker(rand(256));
    %
    % - Then click "Add ROI" to draw an ellipse 
    % - Modify ellipse shape and size if desired
    % - Double-click to confirm the ROI
    % - Optionally add more ROIs
    % - Hit "Export and exit" to quit
    % - ROIs will be saved to base workspace


    properties
        fig         % figure containing the gui
        ax          % axis object for the image
        im          % supplied image data
        overlay     % roi position overlay
        hIm         % handle for the image object
        hOverlay    % handle for the overlay object
        rois        % array of roi structures with masks

        addROI_pb   % pushbutton for adding rois
        exit_pb     % pushbutton to export rois and exit
    end % close properties block


    methods
        function obj = RoiMaker(im)
            % This is the constructor method of the RoiMaker class, it gets 
            % called whenever an instance of the class is created.
            obj.rois = struct('footprint',{});
            obj.im = mat2gray(im) * 3;

            obj.fig = figure;

            obj.overlay = zeros(size(obj.im));

            obj.ax = axes('Position', [0, 0.1, 1, 0.9]);
            obj.hIm = imshow(obj.im);
            hold on

            imRed = zeros(size(im,1), size(im,2), 3);
            imRed(:,:,1) = 1;
            obj.hOverlay = imshow(imRed);
            obj.hOverlay.AlphaData = obj.overlay;

            axis equal off

            % we draw pushbuttons and define the callback functions,
            % which will be called whenever the button is clicked
            obj.addROI_pb = uicontrol('Style', 'pushbutton', ...
                'String', 'Add ROI', ...
                'Units', 'normalized', ...
                'Position', [0, 0, 0.3, 0.1], ...
                'Callback', @obj.add_roi);

            obj.exit_pb = uicontrol('Style', 'pushbutton', ...
                'String', 'Export and exit', ...
                'Units', 'normalized', ...
                'Position', [0.7, 0, 0.3, 0.1], ...
                'Callback', @obj.export_and_exit);
        end

        function obj = add_roi(obj, src, events)
            % create a draggable ellipse
            hRoi = imellipse(obj.ax);

            % wait for user to double click
            wait(hRoi);

            if ~isvalid(hRoi)
                % Do not proceed if user failed to double-click
                return
            end
            % retrieve a mask defined by the boundaries of the ellipse
            mask = createMask(hRoi, obj.hIm);
            delete(hRoi) % get rid of it

            % add a new roi to store the mask
            nRois = numel(obj.rois);
            obj.rois(nRois+1) = struct('footprint', sparse(mask));

            % update overlay display
            obj.overlay(mask) = 0.3;
            obj.hOverlay.AlphaData = obj.overlay;
        end % close add_roi

        function obj = export_and_exit(obj, src, events)
            % assign ROIs to the base workspace and close the figure window
            if length(obj.rois)<1
                fprintf('No ROIs selected\n')
            else
                baseVarName='rois';
                fprintf('Assigning %d ROIs to base workspace variable "%s"\n', length(obj.rois), baseVarName)
                assignin('base', baseVarName, obj.rois)
            end

            close(obj.fig)
        end % close export_and_exit

    end % close methods block

end % close classdef
