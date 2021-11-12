function [Imgout] = imageoverlay(background,foreground,alpha,cax)
% this function takes a background image and a foreground image and applies
% an alpha to the foreground to merge both images
% both the background and foreground can be image stacks of the same length

B=background{1};

ogim = size(B);

imsw = figure;
B1 = figure;
F1 = figure;

for f=1:length(foreground)
    
    figure(B1)
    B=background{f};
    colormap(gray)
    B(B==0)=0/0;
    imagesc(B)
    caxis(cax)
    axis off
    view(2)
    daspect([1,1,1])
    F = getframe;
    bkg{f} = frame2im(F);
        
    figure(F1)
    imshow(foreground{f},[])
    daspect([1,1,1])
    F = getframe;
    fg{f} = frame2im(F);
    
    figure(imsw)
    imshow(imresize(bkg{f},ogim),[])
    hold on
    h = imshow(imresize(fg{f},ogim),[]); % Save the handle; we'll need it later
    hold off
    set(h, 'AlphaData', alpha);
    
    F = getframe;
    Imgout{f} = frame2im(F);
end

close(B1)
close(F1)
close(imsw)