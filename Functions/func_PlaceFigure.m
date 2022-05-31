function [] = func_PlaceFigure(ax,Image)
imagesc(ax,Image)
daspect(ax,[1,1,1])
axis(ax, 'tight');