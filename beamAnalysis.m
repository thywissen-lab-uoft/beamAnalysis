function beamAnalysis

%camera pixel size
% For the camera pixel i would recommend writing a text file in the
% relevant folder where the image is saved
% pxsize = 5.2e-6;  % um per pixel
% pxsize = 3.45e-6; % um per pixel

% Load the image
[filename,mydir,~] = uigetfile({'*.tif';'*.tiff'});
Z = imread(fullfile(mydir,filename));
Z = double(Z);
% Z = double(Z(:,:,1));

% Z=imrotate(Z,15);

% Fit a gaussian
[fout,gof,output]=fitRotatedGaussian(Z);

%% Collect fit outputs
X = 1:size(Z,2);
Y = 1:size(Z,1);

[XX,YY]= meshgrid(X,Y);

s1 = fout.s1;
s2 = fout.s2;
theta = fout.theta;

w1 = 2*s1;
w2 = 2*s2;

A = fout.A;
Xc = fout.Xc;
Yc = fout.Yc;

nbg = fout.nbg;

Zfit = feval(fout,XX,YY);

% 1/e^2 ellipse
tt=linspace(0,2*pi,1000);

Xe = w1*cos(theta)*cos(tt)-w2*sin(theta)*sin(tt);
Ye = w1*sin(theta)*cos(tt)+w2*cos(theta)*sin(tt);

% Major axis
X_a = w1*cos(theta);
Y_a = w1*sin(theta);


% Minor axis
X_b = w2*cos(theta-pi/2);
Y_b = w2*sin(theta-pi/2);

%% Show the output
hF = figure;
hF.Color='w';
hF.Position=[50 50 1800 400];

% Image Plot
ax1 = subplot(131);
h1 = imagesc(X,Y,Z);
% set(gca,'YDir','normal')
axis equal tight
hold on
caxis([0 fout.A*1.2])
colorbar
title(filename,'interpreter','none');

%1/e^2 radius
plot(Xe+Xc,Ye+Yc,'r-')
plot(Xe+Xc,Ye+Yc,'r-')


% major/minor axes
plot(Xc+[-1 1]*X_a,Yc+[-1 1]*Y_a,'r-')
plot(Xc+[-1 1]*X_b,Yc+[-1 1]*Y_b,'r-')

% Residue Plot
ax2 = subplot(132);
h2 = imagesc(X,Y,Z-Zfit);
% set(gca,'YDir','normal')
axis equal tight
colorbar
caxis([-2 2]*nbg);
title('residue');
text(10,10,'look for issues of background! (change fit lims)',...
    'HorizontalAlignment','left','VerticalAlignment','bottom')

% Residue Plot
ax3 = subplot(133);
drawnow;

t = uitable('Parent',hF,'units','normalized','position',ax3.Position);
delete(ax3);

set(t,'RowName',{},'ColumnName',{},'ColumnWidth',{80, 200})

t.Data{1,1}='w1 (px)';
t.Data{2,1}='w2 (px)';
t.Data{3,1}='theta (deg)';
t.Data{4,1}='bg';
t.Data{5,1}='amp';
t.Data{6,1}='Xc';
t.Data{7,1}='Yc';

t.Data{1,2}=w1;
t.Data{2,2}=w2;
t.Data{3,2}=theta*180/pi;
t.Data{4,2}=nbg;
t.Data{5,2}=fout.A;
t.Data{6,2}=Xc;
t.Data{7,2}=Yc;




% 
% for kk=1:length(imgs)
%    Z = double(imgs{kk});
%    X = (1:size(Z,2));
%    Y = (1:size(Z,1));
%     fout=fit2DGauss(X,Y,Z);
%     [xx,yy]=meshgrid(X,Y);    
% 
%     figure
%     subplot(221)
%     imagesc(Z)
%     axis equal tight
% 
%     subplot(222)
%     imagesc(feval(fout,xx,yy));
%     axis equal tight
% 
%     subplot(223)
%     imagesc(Z-feval(fout,xx,yy))
%     axis equal tight
% 
%     w1(kk)=2*fout.s1;
%     w2(kk)=2*fout.s2;
%     theta(kk)=fout.theta;
% 
%     xc(kk)=fout.Xc;
%     yc(kk)=fout.Yc;
% 
%     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%');    
%     disp(filenames{kk});
%     disp(fout);
%     disp(['waist 1 (px) : ' num2str(w1(kk))]);
%     disp(['waist 1 (px) : ' num2str(w2(kk))]);
%     disp(['theta  (rad) : ' num2str(theta(kk))]);
%     disp(['xc (px) : ' num2str(xc(kk))]);
%     disp(['yc (px) : ' num2str(yc(kk))]);
% end


% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%');    

end


