%% Design parameters
disp('Select initial estimation error covariance')
% states: xa, ya, theta_a, xt, yt, theta_t, ut % covariance of initial estimate
xyun0 = input('Enter initial position estimation uncertainty (m): ');
thun0 = input('Enter initial heading estimation uncertainty (rad): ');
spun0 = input('Enter initial speed estimation uncertainty (m/s): ');
px0  = [xyun0 xyun0 thun0 xyun0 xyun0 thun0 spun0].^2/4;
Px0 = diag(px0);

disp('Select process noise covariance')
% process noise covariance
xyunv = input('Enter agent position modeling uncertainty (m/s): ');
thunv = input('Enter agent heading modeling uncertainty (rad/s): ');
yawrv = input('Enter target yaw rate range (rad/s): ');
spunv = input('Enter target speed range (m/s): ');
rv  = [xyunv xyunv thunv xyunv xyunv yawrv spunv].^2/4;
Rv = diag(rv);

disp('Select measurement noise covariance')
% measurement noise covariance: xa, ya, theta_a, xt, yt
xyunn = input('Enter position measurement uncertainty (m): ');
thunn = input('Enter heading measurement uncertainty (rad): ');
spunn = input('Enter speed measurement uncertainty (m/s): ');
rn    = [xyunn xyunn thunn xyunn xyunn].^2/4;
Rn    = diag(rn);

% kinematic constraint
umax = input('Enter max target speed for kinematic constraint (m/s): ');