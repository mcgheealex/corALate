function [] = func_PlaceFigure(ax,Image)
imagesc(ax,Image)
colormap(ax,gray)
daspect(ax,[1,1,1])
axis(ax, 'tight');