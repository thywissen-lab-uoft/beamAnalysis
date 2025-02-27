
function [fout,gof,output]=fitRotatedGaussian(X,Y,Z)

    if nargin==1
        Z = X;
        X = 1:size(Z,2);
        Y = 1:size(Z,1);
    end
    Z=double(Z);

    [s1,s2,theta,Xc,Yc,C] = rotateGaussGuess(Z);

    % Background guess
    nbg=min(min(Z));

    % Amplitdue
    A=max(max(imgaussfilt(Z,.5)))-nbg;

    % Make a mesh grid for fitting
    [xx,yy]=meshgrid(X,Y);    

    % Copy the data
    Z2=Z;xx2=xx;yy2=yy;

    
    % Remove data below threshold
    binds = Z<(0.05*A+nbg);
    Z2(binds)=[];
    xx2(binds)=[];
    yy2(binds)=[];

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
    opt.Upper=[A*1.5 Xc+100 Yc+100 s1*10 s2*10 theta+15*(pi/180) nbg+100];
    opt.Lower=[A*0.5 Xc-100 Yc-100 s1*.2 s2*.2 theta-15*(pi/180) -2];
%opt.Robust='bisquare';
    % Guess Image
    zg=gaussrot(A,Xc,Yc,s1,s2,theta,nbg,xx,yy);

    % Display initial guess  
    gStr=[' guess (Xc,Yc,s1,s2,theta,A,bg)=(' num2str(round(Xc)) ',' ...
        num2str(round(Yc)) ',' num2str(round(s1)) ',' num2str(round(s2)) ',' num2str(round(theta*180/pi,1)) ',' ...
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