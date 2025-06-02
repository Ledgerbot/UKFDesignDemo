%% Load experimental data set
load('USCGA_scenario_1.mat')
Na = size(ukfData,2);

%% Design parameters
% covariance of initial estimate
% states: xa, ya, theta_a, xt, yt, theta_t, ut 
xyun0 = [0.25 0.5 1];  % initial position estimation uncertainty (m)
thun0 = [0.1 0.2 0.4]; % initial heading estimation uncertainty (rad)
spun0 = [0.5 1 2]; % initial speed estimation uncertainty (m/s)

disp('Initial estimation error covariance options: 1 - small, 2 - medium, 3- large')
Px0_choice = input('Select initial estimation error covariance option (1, 2, or 3): ');
if (Px0_choice~=1)&&(Px0_choice~=2)&&(Px0_choice~=3)
    Px0_choice = 2;
end
px0  = [xyun0(Px0_choice) xyun0(Px0_choice) thun0(Px0_choice) xyun0(Px0_choice) ...
    xyun0(Px0_choice) thun0(Px0_choice) spun0(Px0_choice)].^2/4;
Px0 = diag(px0);

disp('Select process noise covariance')
% process noise covariance
% states: xa, ya, theta_a, xt, yt, theta_t, ut 
xyunv = [0.5 1 2];     % agent position modeling uncertainty (m/s)
thunv = [0.2 0.4 0.8]; % agent heading modeling uncertainty (rad/s)
yawrv = [0.4 0.8 1.6]; % target yaw rate range (rad/s)
spunv = [0.5 1 2];     % target speed range (m/s)
disp('Process noise covariance (model uncertainty) options: 1 - small, 2 - medium, 3- large')
Rv_choice = input('Select process noise covariance option (1, 2, or 3): ');
if (Rv_choice~=1)&&(Rv_choice~=2)&&(Rv_choice~=3)
    Rv_choice = 2;
end
rv  = [xyunv(Rv_choice) xyunv(Rv_choice) thunv(Rv_choice) xyunv(Rv_choice) ...
    xyunv(Rv_choice) yawrv(Rv_choice) spunv(Rv_choice)].^2/4;
Rv = diag(rv);

% measurement noise covariance: xa, ya, theta_a, xt, yt
xyunn = [0.25 0.5 1];  % position measurement uncertainty (m)
thunn = [0.1 0.2 0.4]; % heading measurement uncertainty (rad)
spunn = [0.5 1 2];     % speed measurement uncertainty (m/s)
disp('Measurement noise covariance options: 1 - small, 2 - medium, 3- large')
Rn_choice = input('Select process noise covariance option (1, 2, or 3): ');
if (Rn_choice~=1)&&(Rn_choice~=2)&&(Rn_choice~=3)
    Rn_choice = 2;
end
rn    = [xyunn(Rn_choice) xyunn(Rn_choice) thunn(Rn_choice) ...
    xyunn(Rn_choice) xyunn(Rn_choice)].^2/4;
Rn    = diag(rn);

% kinematic constraint
umaxs = [0.25 0.5 1];
disp('Max speed (kinematic constraint) options: 1 - slow, 2 - moderate, 3- fast')
kc_choice = input('Select kinematic constraint option (1, 2, or 3): ');
if (kc_choice~=1)&&(kc_choice~=2)&&(kc_choice~=3)
    kc_choice = 2;
end
umax = umaxs(kc_choice);
%% Initial conditions
xh0 = zeros(7,Na);
% assume initial speed and heading are zero
for i = 1:Na
    xh0(1,i)   = Xa(1,i);
    xh0(2,i)   = Ya(1,i);
    xh0(3,i)   = Yaw(1,i);
    xh0(4:5,i) = loc_tar0(:,i);
end

%% likelihood of broadcasting
threshold = 1;
broadcast = rand(Na,Nsteps);
B = broadcast<=threshold;
NB = broadcast > threshold;
tmp = broadcast;
tmp(B) = 1;
tmp(NB) = 0;
broadcast = tmp;
%% UKF
agents = create3_ukf_exp_fusion_kincon(xh0,Px0,Rv,Rn,Xa,Ya,Yaw,YawRt,U,Y,Ts,Na,umax,broadcast);
%% Plots - scenario
agnts_view = sum(~isnan(Y),2);
trans = find(diff(agnts_view )~=0);

if length(agents(1,1).xh_k(4,:))>1000
    indx = 1:20:length(agents(1,1).xh_k(4,:));
elseif length(agents(1,1).xh_k(4,:))>100
    indx = 1:2:length(agents(1,1).xh_k(4,:));
else
    indx = 1:length(agents(1,1).xh_k(4,:));
end

figure
ax1(1) = subplot(311);
hold on
plot(time(indx),Xt(indx,1)-agents(1,1).xh_k(4,indx)','r','LineWidth',2)
errorbar(time(indx),Xt(indx,1)-agents(1,1).xh_k(4,indx)',2*sqrt(squeeze(agents(1,1).Px_k(4,4,indx))),'LineStyle','none','Color','b'),grid %,'LineWidth',1
hold off
ylabel('Error X (m)')
ax1(2) = subplot(312);
hold on
plot(time(indx),Yt(indx,1)-agents(1,1).xh_k(5,indx)','r','LineWidth',2)
errorbar(time(indx),Yt(indx,1)-agents(1,1).xh_k(5,indx)',2*sqrt(squeeze(agents(1,1).Px_k(5,5,indx))),'LineStyle','none','Color','b'),grid
hold off
ylabel('Error Y (m)')
ax1(3) = subplot(313);
plot(time,ra,'s',time(1),0,time(end),0),grid
ylim([0 3])
ylabel('Rel. agent angle (rad)')
xlabel('Time (s)')
linkaxes(ax1,'x')

figure
ax2(1) = subplot(211);
hold on
grid
plot(time(indx),Xt(indx,1)-agents(1,1).xh_k(4,indx)','LineWidth',2)
plot(time(indx),Yt(indx,1)-agents(1,1).xh_k(5,indx)','LineWidth',2)
ylabel('Error (m)')
hold off
legend('X','Y','Location','best')
ax2(2) = subplot(212);
plot(time,ra,'s',time(1),0,time(end),0),grid
ylim([0 3])
ylabel('Rel. agent angle (rad)')
xlabel('Time (s)')
linkaxes(ax2,'x')