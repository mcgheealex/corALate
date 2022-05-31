function [] = func_PlaceVol(ax,Volume,zscale,viewType)
reset(ax)
axis(ax,'off')

figure('visible','off');
[x,y,z] = size(Volume);
[X,Y,Z] = meshgrid(1:1:x,1:1:y,1:1:z);
V2 = medfilt3(Volume);

thresh = max(Volume,[],'all')*0.1;
V2(V2<thresh) = 0;

X(V2==0)=[];
Y(V2==0)=[];
Z(V2==0)=[];

projection_xy = func_project(Volume,'XY');
projection_xz = func_project(Volume,'XZ');
projection_yz = func_project(Volume,'YZ');

figure('visible','off');
subtightplot(2,2,1,0,0,0)
imagesc(projection_xz)
axis off
colormap('bone')
title('xz')
subtightplot(2,2,2,0,0,0)
imagesc(projection_xy)
axis off
colormap('bone')
title('xy')
subtightplot(2,2,4,0,0,0)
imagesc(projection_yz)
axis off
colormap('bone')
title('yz')
subtightplot(2,2,3,0,0,0)
plot3(X,Y,Z*zscale,'.r')
axis off
title('vol')

f = getframe(gcf);
im2D = f.cdata;

switch viewType
    case '3D'
        plot3(ax,X,Y,Z*zscale,'.r')
        axis(ax,'off')
    case 'Multi2D'
        imshow(im2D,[],'Parent',ax);
end

close all

end
