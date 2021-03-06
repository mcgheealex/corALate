function [Imoverlay,ImRaw] = PlotstrainQuadtree2(U,F,coordinatesFEM,elementsFEM,CurrentImg,DICpara,CAXIS)
%PLOTSTRAINQUADTREE: to plot DIC solved strain fields on a quadtree mesh 
% and overlaid with the original DIC images 
%   [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
%    strain_maxshear,strain_vonMises] = PlotstrainQuadtree(U,F,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
% ----------------------------------------------
%
%   INPUT: F                    DIC solved deformation gradient tensor
%          coordinatesFE        FE mesh coordinates
%          elementsFEM          FE mesh elements
%
%   OUTPUT: strain_exx              strain xx-compoent
%           strain_exy              strain xy-compoent
%           strain_eyy              strain yy-compoent
%           strain_principal_max    max principal strain on the xy-plane
%           strain_principal_min    min principal strain on the xy-plane
%           strain_maxshear         max shear strain on the xy-plane
%           strain_vonMises         equivalent von Mises strain
%
%   Plots:       
%       1) strain sxx
%       2) strain sxy
%       3) strain syy
%       4) max principal strain on the xy-plane 
%       5) min principal strain on the xy-plane
%       6) max shear strain on the xy-plane
%       7) equivalent von Mises strain
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImRaw = [];

%% Initialization
% warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
um2px = DICpara.um2px; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

disp_u = U(1:2:end); disp_v = U(2:2:end);
coordinatesFEMWorldDef = [coordinatesFEM(:,1)+Image2PlotResults*disp_u, coordinatesFEM(:,2)+Image2PlotResults*disp_v];


%% Compute strain components

u_x = F(1:4:end); v_x = F(2:4:end);
u_y = F(3:4:end); v_y = F(4:4:end);

strain_exx = u_x; 
strain_exy = 0.5*(v_x+u_y);
strain_eyy = v_y;

strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
% Principal strain
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
% equivalent von Mises strain
strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain exx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure('visible','off'); ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_exx,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo); caxis([0,0.5]) % D Sample 
% colormap(turbo); caxis([-0.1,0.02]) % foam
% colormap(turbo); caxis([-0.004,0]); % Sample 12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(CAXIS(1,:))
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 frame = getframe(fig1); 
im = frame2im(frame);
Imoverlay.E_xx = im;
ImRaw.E_xx = reshape(strain_exx,[20,520/20]);
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) Strain exy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure('visible','off'); ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_exy,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo); caxis([-0.08,0.08]) % D Sample 
% colormap(turbo); caxis([-0.06,0.06]) % foam
% colormap(turbo); caxis([-0.008,0.008]); % Sample 12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(CAXIS(2,:))
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 frame = getframe(fig1); 
im = frame2im(frame);
Imoverlay.E_xy = im;
ImRaw.E_xy = reshape(strain_exy,[20,520/20]);
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain eyy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure('visible','off'); ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_eyy,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo); caxis([-0.15,0]) % D Sample 
% colormap(turbo); caxis([-0.05,0.2]) % foam
% colormap(turbo); caxis([-0.002,0.017]); % Sample 12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(CAXIS(3,:))
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 frame = getframe(fig1); 
im = frame2im(frame);
Imoverlay.E_yy = im;
ImRaw.E_yy = reshape(strain_eyy,[20,520/20]);
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 4) Strain e_principal_max ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure('visible','off'); ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_principal_max,'NoEdgeColor');

max(strain_principal_max,[],"all")

set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo);  caxis auto; % D Sample 
% colormap(turbo); caxis auto % foam
% colormap(turbo); caxis([0,0.02]); % Sample 12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(CAXIS(4,:))
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

frame = getframe(fig1); 
im = frame2im(frame);
Imoverlay.E_pmax = im;
ImRaw.E_pmax = reshape(strain_principal_max,[20,520/20]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 5) Strain e_principal_min ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure('visible','off'); ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_principal_min,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo);  caxis auto; % D Sample 
% colormap(turbo); caxis auto % foam
% colormap(turbo); caxis([-0.008,0]); % Sample 12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(ax2,CAXIS(5,:))
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

frame = getframe(fig1); 
im = frame2im(frame);
Imoverlay.E_pmin = im;
ImRaw.E_pmin = reshape(strain_principal_min,[20,520/20]); 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 6) Strain e_max_shear ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure('visible','off'); ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_maxshear,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo);  caxis auto; % D Sample 
% colormap(turbo); caxis auto % foam
% colormap(turbo); caxis([0,0.011]); % Sample 12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(ax2,CAXIS(6,:))
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

frame = getframe(fig1); 
im = frame2im(frame);
Imoverlay.E_ms = im;
ImRaw.E_ms = reshape(strain_maxshear,[20,520/20]); 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 7) von Mises equivalent strain ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure('visible','off'); ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_vonMises,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(turbo);  caxis auto; % D Sample 
% colormap(turbo); caxis auto % foam
% colormap(turbo); caxis([0,0.025]); % Sample 12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(ax2,CAXIS(7,:))
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

frame = getframe(fig1); 
im = frame2im(frame);
Imoverlay.E_vm = im;
ImRaw.E_vm = reshape(strain_vonMises,[20,520/20]);

end
 
 
