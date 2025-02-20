function  ALS_propagate
disp(' ');

%% Functions
lambda = 1064e-9;

q2w = @(q) sqrt(-lambda/(pi*imag(1/q)));
q2zr = @(q) imag(q);
q2w0 = @(q) sqrt(q2zr(q)*lambda/pi);
q2z0 = @(q) -real(q);


%% Output of laser head
w0 = 0.45e-3;
z0  = 1.4;
q0 = (0 - z0) + 1i*pi*w0^2/lambda;

dispQ(q0);

%% Travel to first telescope
d = 525e-3;
q1 = q0 + d;

dispQ(q1);


%% First Telescope
f1 = 150e-3;
f2 = -75e-3;
d = 63e-3;

% f1 = 100e-3;
% f2 = -75e-3;
% d = 19e-3;

% Lens, distance, elns
q2a = q1/(-q1/f1+1);
q2b = q2a + d;
q2c = q2b/(-q2b/f2+1);

% Output of telescope
q2 = q2c;

dispQ(q2);

%% Travel through AOM and to second telescope
% d = 500e-3;
% q3 = q2 + d;
% 
% dispQ(q3);


%%

function dispQ(q)

w = q2w(q);
zr = q2zr(q);
w00 = q2w0(q);
z00 = q2z0(q);

disp(['D(um) : ' num2str(2*w*1e6) '; D0(um) : ' num2str(2*w00*1e6) ...
    '; zr(mm) : ' num2str(zr*1e3) '; z0(mm) ' num2str(z00*1e3)]);

end

end

