 w1=5.2*[  97.5306;
   44.0824;
   40.0840;
   43.1724;
   81.9071;
   40.0840];

w2=5.2*[   91.7698;
   41.5608;
   38.4898;
   39.8024;
   79.6997;
   38.4898];

wbar = sqrt(w1.*w2);

z = [14 480.6 389.65 357.9 142 389.65]';

wbar=wbar*1e-6;
z=z*1e-3;

figure
plot(z,wbar,'ko');
xlabel('mm')
ylabel('f pump waist (um)');
lambda = 767e-9;

myfit = fittype(@(w0,z0,z) w0*sqrt(1+((z-z0)/(pi*w0^2/lambda)).^2),...
    'independent',{'z'},'coefficients',{'w0','z0'});
fitopt = fitoptions(myfit);
fitopt.StartPoint = [200e-6 .400];
fout=fit(z,wbar,myfit,fitopt);
hold on

zz=linspace(min(z),max(z),1e3);

plot(zz,feval(fout,zz),'r-');

disp(fout);
