% Section 3* for quadtree mesh

%%
% ====== Generate an image mask in the reference image ======
ImgRef = (imread(file_name{1})); 
try ImgRef = rgb2gray(ImgRef); catch; end % Change the image to grayscale image if it is rgb image

%% ====== Find center of the hole ======
figure, imshow(ImgRef);  % Read reference image
title('Click more than 5 points along the hole boundary, then press -Enter- key','fontweight','normal');
fprintf('--- Define hole geometry ---  \n');
fprintf('Click more than 5 points along the hole boundary,  \n');
fprintf('then press -Enter- key. \n');
[gridx,gridy] = ginput(); ParHole = CircleFitByTaubin([gridx,gridy]); % ParHole = [xCenCoord, yCenCoord, rad]; ParHole=[194.146,542.51,52.072];
fprintf('--- Define hole geometry: Done! ---  \n');


%% ====== Get the grayscale value threshold and construct the image mask ======
ImgRefGaussFilt = imgaussfilt(imgaussfilt(ImgRef,0.5),0.5); % Gaussian filter

close all; figure, imshow(ImgRefGaussFilt);  % Read reference image
title('Click more than 3 points inside the hole, then press -Enter- key','fontweight','normal');
fprintf('--- Define image mask file ---  \n');
fprintf('Click more than 3 points inside the hole \n');
fprintf('Then press -Enter- key. \n');
[gridx,gridy] = ginput();  ImgRefMaskThresholdInside = 0;
for tempi = 1:length(gridx)
    ImgRefMaskThresholdInside = ImgRefMaskThresholdInside + double(ImgRefGaussFilt(round(gridy(tempi)),round(gridx(tempi))));
end
ImgRefMaskThresholdInside = ImgRefMaskThresholdInside/length(gridx);
% For example, ImgRefMaskThreshold = 20; in the imgFiles_Quadtree_Sample12 exampple

close all; figure, imshow(ImgRefGaussFilt);  % Read reference image
title('Click more than 3 points outside the hole, then press -Enter- key','fontweight','normal');
fprintf('--- Define image mask file ---  \n');
fprintf('Click more than 3 points outside the hole \n');
fprintf('Then press -Enter- key. \n');
[gridx,gridy] = ginput();  ImgRefMaskThresholdOutside = 0;
for tempi = 1:length(gridx)
    ImgRefMaskThresholdOutside = ImgRefMaskThresholdOutside + double(ImgRefGaussFilt(round(gridy(tempi)),round(gridx(tempi))));
end
ImgRefMaskThresholdOutside = ImgRefMaskThresholdOutside/length(gridx);
% For example, ImgRefMaskThreshold = 20; in the imgFiles_Quadtree_Sample12 exampple

fprintf('\n'); fprintf('Is the hole darker or brighter than the sample surface?  \n');
fprintf('    0: Hole is darker; \n');
fprintf('    1: Hole is brighter; \n');
prompt = 'Input here: '; HoleDarkerOrBrighter = input(prompt);

if HoleDarkerOrBrighter == 0
    ImgRefMask = logical(ImgRefGaussFilt > 0.5*(ImgRefMaskThresholdInside+ImgRefMaskThresholdOutside));
elseif HoleDarkerOrBrighter == 1
    ImgRefMask = logical(ImgRefGaussFilt < 0.5*(ImgRefMaskThresholdInside+ImgRefMaskThresholdOutside));
end

[tempxx,tempyy] = ndgrid(1:size(ImgRefMask,1), 1:size(ImgRefMask,2));
dist2HoleCenter = sqrt( (tempxx-ParHole(2)).^2 + (tempyy-ParHole(1)).^2 );
[row,col] = find(dist2HoleCenter < ParHole(3));
for tempi = 1:length(row)
   ImgRefMask(row(tempi),col(tempi)) = 0; 
end
figure, imshow(ImgRefMask); title('Image Mask')

removeobjradius =  round(0.2*ParHole(3)); % remove all object containing fewer than hole Rad
se = strel('disk',removeobjradius); ImgRefMask = imclose(ImgRefMask,se); imshow(ImgRefMask);


%% Deal with outside  DICpara.gridxyROIRange
gridx = DICpara.gridxyROIRange.gridx;
gridy = DICpara.gridxyROIRange.gridy;

ImgRefMask(:, 1:gridx(1)) = 1;
ImgRefMask(:, gridx(2):end) = 1;
ImgRefMask(1:gridy(1), :) = 1;
ImgRefMask(gridy(2):end, :) = 1;
figure, imshow(ImgRefMask);
 

%% ====== Store ImgRefMask ======
DICpara.ImgRefMask = ImgRefMask;


