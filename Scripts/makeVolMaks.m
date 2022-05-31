activebutton = app.ButtonGroup_ROI_3D.SelectedObject.Text;

% use the current image to get the size and to use it as a
% background
Im = ones([1000,1000]); % all 3D vol masks strat as 1000,1000 size
[m,n,k] = size(app.ROIMask); % get the size of the volume

% initialize the masks
add2image = zeros(size(Im));
subfromimage = zeros(size(Im));

% go through the current axis and see if there are any shapes
% drawn on it, if so use those shapes to make a mask
for i=1:numel(app.UIAxes_ROI_3D.Children)
    % the children are the rectangles, circles etc select one
    % of them
    ax = app.UIAxes_ROI_3D.Children(i);
    % need to use a try and catch because some elements of the
    % axis may not have labels
    try
        if isequal(ax.Label,'add') % if one of the objects has this label
            mask = createMask(ax,Im);
            add2image = add2image + mask;
            add2image(add2image>1)=1; % if a mask overlaps then we dont care
        elseif isequal(ax.Label,'subtract')
            mask = createMask(ax,Im);
            subfromimage = subfromimage + mask;
            subfromimage(subfromimage>1)=1; % if a mask overlaps then we dont care
        end
    catch
        % do nothing
    end
end


switch activebutton
    
    case 'XY projection'
        % resize the add / subtract image to the roi size
        add2image = imresize(add2image,[m,n]);
        subfromimage = imresize(subfromimage,[m,n]);
        % multiply one projected axis by the add/subtract
        
        for kk=1:k
            % from all of the add and subtract rectangles, make a mask
            app.ROIMask(:,:,kk) = app.ROIMask(:,:,kk) + add2image;
            app.ROIMask(:,:,kk) = app.ROIMask(:,:,kk) - subfromimage;
            
        end
        app.ROIMask(app.ROIMask>0.5) = 1;
        app.ROIMask(app.ROIMask<0.5) = 0;
        % set the ROI axes to the first image
        IM = func_project(app.AllImages{1},'XY');
        ROI = func_project(app.ROIMask,'XY');
        doplot = true;
    case 'YZ Projection'
        % resize the add / subtract image to the roi size
        add2image = imresize(add2image,[n,k]);
        subfromimage = imresize(subfromimage,[n,k]);
        % multiply one projected axis by the add/subtract
        a2i(1,:,:) = add2image;
        sfi(1,:,:) = subfromimage;
        for mm=1:m
            % from all of the add and subtract rectangles, make a mask
            app.ROIMask(mm,:,:) = app.ROIMask(mm,:,:) + a2i;
            app.ROIMask(mm,:,:) = app.ROIMask(mm,:,:) - sfi;
        end
        app.ROIMask(app.ROIMask>0.5) = 1;
        app.ROIMask(app.ROIMask<0.5) = 0;
        % set the ROI axes to the first image
        IM = func_project(app.AllImages{1},'YZ');
        ROI = func_project(app.ROIMask,'YZ');
        doplot = true;
    case 'XZ projection'
        % resize the add / subtract image to the roi size
        add2image = imresize(add2image,[m,k]);
        subfromimage = imresize(subfromimage,[m,k]);
        a2i(:,1,:) = add2image;
        sfi(:,1,:) = subfromimage;
        % multiply one projected axis by the add/subtract
        for nn=1:n
            % from all of the add and subtract rectangles, make a mask
            app.ROIMask(:,nn,:) = app.ROIMask(:,nn,:) + a2i;
            app.ROIMask(:,nn,:) = app.ROIMask(:,nn,:) - sfi;
        end
        app.ROIMask(app.ROIMask>0.5) = 1;
        app.ROIMask(app.ROIMask<0.5) = 0;
        % set the ROI axes to the first image
        IM = func_project(app.AllImages{1},'XZ');
        ROI = func_project(app.ROIMask,'XZ');
        doplot = true;
    case 'Volume'
        % do nothing
        doplot = false;
end

if doplot
    ROIfig = imageoverlay_Vol(IM,ROI,0.5);
    func_PlaceFigure(app.UIAxes_ROI_3D,ROIfig);
end