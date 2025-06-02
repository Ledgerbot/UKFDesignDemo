function f = create3_exp(z,u,w,T)
% Includes kinematic constraints on target speed and yaw rate
%% States: z = x, y, theta, L, phi 
% xa      = z(1);
% ya      = z(2);
theta_a = z(3);
% xt      = z(4);
% yt      = z(5);
theta_t = z(6);
ut      = z(7);

%% Inputs: r (known yaw rate), w (process noise)
% Agent (known)
ra  = u(1); % yaw rate
ua  = u(2); % forward speed

% Target
rt  = w(1); % yaw rate process noise
utw = w(2); % speed estimation process noise

%% CT model
% Kinematic constraints
u_max = 0.5; % m/s
r_max = 3; % rad/s

% Agent
if ua > u_max
    ua = u_max;
elseif ua < -u_max
    ua = -u_max;
end

if ra > r_max
    ra = r_max;
elseif ra < -r_max
    ra = -r_max;
end

% Target
if ut > u_max
    ut = u_max;
elseif ut < -u_max
    ut = -u_max;
end

if rt > r_max
    rt = r_max;
elseif rt < -r_max
    rt = -r_max;
end


fc    = zeros(numel(z),1);
fc(1) = ua*cos(theta_a);
fc(2) = ua*sin(theta_a);
fc(3) = ra;
fc(4) = ut*cos(theta_t);
fc(5) = ut*sin(theta_t);
fc(6) = rt;
fc(7) = utw;

%% DT model
f = z + fc*T;
end