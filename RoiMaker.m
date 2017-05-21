classdef RoiMaker < handle
    % A simple GUI to interactively place elliptical ROIs on an image
    properties
        fig         % figure containing the gui
        ax          % axis object for the image
        im          % supplied image data
        overlay     % roi position overlay
        him         % handle for the image object
        hoverlay    % handle for the overlay object
        rois        % array of roi structures with masks
        
        addroipb    % pushbutton for adding rois
        exitpb      % pushbutton to export rois and exit
    end
    
    methods
        function obj = RoiMaker(im)
            % this is the constructor method of the RoiMaker class
            % it gets called whenever an instance of the class is created
            obj.rois = struct('footprint',{});
            obj.im = mat2gray(im) * 3;

            obj.fig = figure;
            
            obj.overlay = zeros(size(obj.im));
            
            obj.ax = axes('Position', [0 .1 1 .9]);
            obj.him = imshow(obj.im); hold on;
            
            imred = zeros(size(im,1), size(im,2), 3);
            imred(:,:,1) = 1;
            obj.hoverlay = imshow(imred);
            obj.hoverlay.AlphaData = obj.overlay;  
            
            axis equal off;
            
            % we draw pushbuttons and define the callback functions,
            % which will be called whenever the button is clicked
            obj.addroipb = uicontrol('Style', 'pushbutton', ...
                'String', 'Add ROI',...
                'Units', 'normalized',...
                'Position', [0 0 .3 .1],...
                'Callback', @obj.add_roi);
            
            obj.exitpb = uicontrol('Style', 'pushbutton', ...
                'String', 'Export and exit',...
                'Units', 'normalized',...
                'Position', [0.7 0 .3 .1],...
                'Callback', @obj.export_and_exit);
        end
        
        function obj = add_roi(obj, src, events)
            % create a draggable ellipse
            hroi = imellipse(obj.ax);
            % wait for user to double click
            wait(hroi);
            % retrieve a mask defined by the boundaries of the ellipse
            mask = createMask(hroi, obj.him);
            delete(hroi);   % get rid of it
            % add a new roi to store the mask
            nRois = numel(obj.rois);
            obj.rois(nRois+1) = struct('footprint', sparse(mask));
            
            % update overlay display
            obj.overlay(mask) = .3;
            obj.hoverlay.AlphaData = obj.overlay;
        end
        
        function obj = export_and_exit(obj, src, events)
            % assign rois to the base workspace and close
            assignin('base', 'rois', obj.rois);
            close(obj.fig);
        end
    end
end
            
            
