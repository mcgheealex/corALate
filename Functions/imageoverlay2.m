function [imgout] = imageoverlay2(background,foreground,alpha)
% this function takes a background image and a foreground image and applies
% an alpha to the foreground to merge both images
% both the background and foreground can be image stacks of the same length

ogim = size(background,1,2);
% AxesH = axes('Units', 'pixels', 'position', [10, 10, ogim(1), ogim(2)], ...
%              'Visible', 'off');
         
B1 = figure('visible','off');
B1.Position = [1,1,ogim(2)+100,ogim(1)+100];
colormap(gray)
background(background==0)=0/0;
imshow(background)
axis off
view(2)
daspect([1,1,1])
F = getframe;
bkg = rgb2gray(frame2im(F));

F1 = figure('visible','off');
F1.Position = [1,1,ogim(2)+100,ogim(1)+100];
imshow(foreground,[])
daspect([1,1,1])
F = getframe;
fg = rgb2gray(frame2im(F));

imsw = figure('visible','off');
imsw.Position = [1,1,ogim(2)+100,ogim(1)+100];
imshow(imresize(bkg,ogim(:,:,1)),[])
hold on
h = imshow(imresize(fg,ogim),[]); % Save the handle; we'll need it later
hold off
set(h, 'AlphaData', alpha);

F = getframe;
imgout = imresize(frame2im(F),ogim);

close(B1)
close(F1)
close(imsw)