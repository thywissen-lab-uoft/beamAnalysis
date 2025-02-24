function [s1,s2,theta,muX,muY,C] = rotateGaussGuess(img)
% Author : CJ Fujiwara
%
% Make an initial guess for fitting a rotated gaussian.  This code looks at
% the eigenvectors of the covariance matrix.

%% Pre-pocessing
img = double(img);              % Make it double precision
bg = min(min(double(img)));     % Find background
image_bgremoved = img - bg;     % Remove background

% Normalize image
image_normalized = image_bgremoved/sum(image_bgremoved,'all');

% Compute X and Y vectors
X = 1:size(image_normalized,2);
Y = 1:size(image_normalized,1);
[XX,YY]=meshgrid(X,Y);

%% Compute Center of Mass
% Find Center of mass
muX = sum(XX.*image_normalized,'all');
muY = sum(YY.*image_normalized,'all');
%% Compute Covariances
% Variance
varXX = sum((XX-muX).^2.*image_normalized,'all');
varYY = sum((YY-muY).^2.*image_normalized,'all');

% Covariance
varXY = sum((XX-muX).*(YY-muY).*image_normalized,'all');

% Covariance matrix
C = [varXX varXY;
    varXY varYY];

%% Solve Eigenvalue

% Eigenvalues of covariance matrix
[U,D] = eig(C);
D=diag(D);

%% Construct output

v1 = U(:,1);    % Small Vector
v2 = U(:,2);    % Big Vector

% Angle of small waist
theta = atan2(v1(2),v1(1));

% Construct waists
s1 = sqrt(D(1)); % Small waist
s2 = sqrt(D(2)); % Big Waist

end

