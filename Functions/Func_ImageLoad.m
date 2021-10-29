%% load images from folder
function [ImageFiles,baseName] = Func_ImageLoad(varargin)
% make a struct with each variable being an option to parse with the
% function parseargs
X.multi = 'on';
X = parseargs(X,varargin{:});

try
    % set the types of images we can see in the folder
    filterspec = {'*.jpg;*.jpeg;*.tif;*.tiff;*.bmp;*.png;*.jp2'};
    % switch between either multiselect on or off
    switch X.multi
        case 'on'
            % open a UI for browsing for images
            [baseName, folder] = uigetfile(filterspec, 'MultiSelect', 'on');
        case 'off'
            [baseName, folder] = uigetfile(filterspec, 'MultiSelect', 'off');
        otherwise
            % should never happen
            disp('Error in func_ImageLoad (no case found)')
    end
    
    % Make sure user didn't cancel uigetfile dialog
    if (ischar(folder))
       fname = fullfile(folder, baseName);
       ImageFiles = fname;
    end
   
    
catch
    fprintf('Error in func_ImageLoad')
    ImageFiles = [];
    baseName = [];
end
