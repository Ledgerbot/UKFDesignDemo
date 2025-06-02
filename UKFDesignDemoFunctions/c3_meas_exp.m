% x = agent states, target (x,y)
% v: state measurement noise
function z = c3_meas_exp(x,v)
xa      = x(1) + v(1);
ya      = x(2) + v(2);
theta_a = x(3) + v(3);
xt      = x(4) + v(4);
yt      = x(5) + v(5);
% theta_t = x(6) + v(6);

dx    = xt - xa;
dy    = yt - ya;
phi   = atan2(dy,dx);

z1 = [xa;ya;theta_a];
z  = [z1;phi-theta_a];

end