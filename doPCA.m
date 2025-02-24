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