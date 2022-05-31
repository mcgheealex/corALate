function [imgout] = imageoverlay_Vol(background,foreground,alpha)

imgout = double(background).*double(foreground*alpha);

end