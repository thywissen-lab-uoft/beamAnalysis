function [w1,w2,theta,xc,yc]=beamAnalysis

ext = 'tif';
newdir=uigetdir(pwd,'select folder of images');

pxsize = 5.2E-6; % um per pixel

if ~newdir
    return;
end

filenames = dir([newdir filesep '*.' ext]);
filenames={filenames.name};
imgs={};
for kk=1:length(filenames)
    A = imread(fullfile(newdir,filenames{kk}));
    A = A(:,:,1);
    imgs{kk}=A;    
end

w1=zeros(length(imgs),1);
w2=zeros(length(imgs),1);
xc=zeros(length(imgs),1);
yc=zeros(length(imgs),1);
theta=zeros(length(imgs),1);

for kk=1:length(imgs)
   Z = double(imgs{kk});
   X = (1:size(Z,2));
   Y = (1:size(Z,1));
    fout=fit2DGauss(X,Y,Z);
    [xx,yy]=meshgrid(X,Y);    


    figure
    subplot(221)
    imagesc(Z)
    axis equal tight

    subplot(222)
    imagesc(feval(fout,xx,yy));
    axis equal tight

    subplot(223)
    imagesc(Z-feval(fout,xx,yy))
    axis equal tight

    w1(kk)=2*fout.s1;
    w2(kk)=2*fout.s2;
    theta(kk)=fout.theta;
    
    xc(kk)=fout.Xc;
    yc(kk)=fout.Yc;
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%');    
    disp(filenames{kk});
    disp(fout);
    disp(['waist 1 (px) : ' num2str(w1(kk))]);
    disp(['waist 1 (px) : ' num2str(w2(kk))]);
    disp(['theta  (rad) : ' num2str(theta(kk))]);
    disp(['xc (px) : ' num2str(xc(kk))]);
    disp(['yc (px) : ' num2str(yc(kk))]);

end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%');    


end

function fout=fit2DGauss(X,Y,Z)

    PCA = doPCA(Z);

% Rotated Gaussian Analysis
    % Very experimental
    % Z is a mxn matrix to which we want to fit an elliptic gaussian            

    % Angle of rotation     
    theta=atan(PCA.PCA(2,1)/PCA.PCA(1,1));

    % Gaussian Radius
    s1=PCA.Radii(1);
    s2=PCA.Radii(2);

    % Center point
    Xc=PCA.Mean(1);
    Yc=PCA.Mean(2);
    % Amplitdue
    A=max(max(imgaussfilt(Z,.5)));

    % Background guess
    nbg=min(min(Z));

    % Make a mesh grid for fitting
    [xx,yy]=meshgrid(X,Y);    

    % Copy the data
    Z2=Z;xx2=xx;yy2=yy;

    % https://en.wikipedia.org/wiki/Gaussian_function
    % But we add a minus sign to make it counter clockwise angle
    % When theta=0 s1 is on the x axis
    gaussrot=@(A,Xc,Yc,s1,s2,theta,nbg,xx,yy) A*exp(-( ...
        (cos(theta)^2/(2*s1^2)+sin(theta)^2/(2*s2^2))*(xx-Xc).^2 + ...
         2*(sin(2*theta)/(4*s1^2) - sin(2*theta)/(4*s2^2))*(xx-Xc).*(yy-Yc) + ...
         (sin(theta)^2/(2*s1^2)+cos(theta)^2/(2*s2^2))*(yy-Yc).^2))+nbg;   

    myfit=fittype(@(A,Xc,Yc,s1,s2,theta,nbg,xx,yy) gaussrot(A,Xc,Yc,s1,s2,theta,nbg,xx,yy),...
        'independent',{'xx','yy'},'coefficients',{'A','Xc','Yc','s1','s2','theta','nbg'});
    opt=fitoptions(myfit);

    opt.StartPoint=[A Xc Yc s1 s2 theta nbg];
    opt.Upper=[A*1.3 Xc+50 Yc+50 s1*5.5 s2*5.5 theta+5*(pi/180) 200];
    opt.Lower=[A*0.7 Xc-50 Yc-50 s1*.5 s2*.5 theta-5*(pi/180) -2];

    % Display initial guess  
    gStr=[' guess (Xc,Yc,s1,s2,theta,A,bg)=(' num2str(round(Xc)) ',' ...
        num2str(round(Yc)) ',' num2str(round(s1)) ',' num2str(round(s2)) ',' num2str(round(theta*pi/180,1)) ',' ...
        num2str(A,'%.e') ',' num2str(round(nbg)) ')' ];     
    fprintf([gStr ' ... ']);

    % Perform the fit
    t1=now;
    [fout,gof,output]=fit([xx2(:) yy2(:)],Z2(:),myfit,opt);
    t2=now;            

    % Fit String
    fStr=['(' num2str(round(fout.Xc)) ',' ...
        num2str(round(fout.Yc)) ',' num2str(round(fout.s1)) ',' num2str(round(fout.s2)) ',' num2str(round(fout.theta*pi/180,1)) ',' ...
        num2str(fout.A,'%.e') ',' num2str(round(fout.nbg)) ')' ];     
    fprintf([' fit ' fStr]);

    % dispaly Summary
    disp([' (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);

     
end

function PCA = doPCA(Z)
    PCA=struct;    
   
    % Smooth and normalize data
    B=imgaussfilt(Z,2);
    B=B/max(max(B));
    
    % Find 30%-40% 1/e=.36
    C=(B>.32).*(B<.4);
    [y,x]=find(C);
    
    pp=[x y];

    
    % Do PCA
    [cc,score,latent,tsquared,explained,mu]=pca(pp);

    % Gaussian radius
    s1=sqrt(latent(1)); % exp(-1/2) guess
    s2=sqrt(latent(2)); % exp(-1/2) guess
    
    % X and Y center points
    xc=mu(1);
    yc=mu(2);
    
    % Normalized principal vector
    v1=cc(1,:);
    v2=cc(2,:);    
    
    PCA.Mean=[xc yc];
    PCA.Radii=[s1 s2];
    PCA.PCA=cc;

end