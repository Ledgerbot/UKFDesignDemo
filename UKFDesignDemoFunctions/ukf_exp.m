function [x,P,x1,P1,P2]=ukf_exp(fstate,x,P,hmeas,z,Q,R,Uk,meas,ro,Ts)
% UKF   Unscented Kalman Filter for nonlinear dynamic systems
% [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: fanction handle for h(x)
%           z: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance
%
% Example:
%{
n=3;      %number of state
q=0.1;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation
s= [1;1;0]; %[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covraiance
N=100;                                     % total dynamic steps
xV = zeros(n,N);          %estmate        % allocate memory
sV = zeros(n,N);          %actual
zV = zeros(1,N);
for k=1:N
  z = h(s) + r*randn;                     % measurments
  sV(:,k)= s;                             % save actual state
  zV(k)  = z;                             % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  s = f(s) + q*randn(3,1);                % update process 
end
for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
end
%}
% Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
% Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
% pp.401-422, 2004. 
%
% By Yi Cao at Cranfield University, 04/01/2008
%
Lx = size(P,1);
Lq = size(Q,1);
Lr = size(R,1);                                 
L = Lx + Lq  + Lr;                     %numer of states + dim of process noise + dim of meas noise
m=numel(z);                                 %numer of measurements
alpha=1e-3;                                 %default, tunable
ki=0;                                       %default, tunable
beta=2;                                     %default, tunable
lambda=alpha^2*(L+ki)-L;                    %scaling factor
c=L+lambda;                                 %scaling factor
Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
c=sqrt(c);
Pf = zeros(L,L);
Pf(1:Lx,1:Lx) = P;
Pf(Lx+1:Lx+Lq,Lx+1:Lx+Lq) = Q;
Pf(Lx+Lq+1:L,Lx+Lq+1:L) = R;
X=sigmas([x;zeros(Lq+Lr,1)],Pf,c);           %sigma points around x
[x1,X1,P1,X2]=ut_knwn_inpt_exp(fstate,X,Wm,Wc,Lx,Q,Uk,Ts);  %unscented transformation of process
% X1=sigmas(x1,P1,c);                         %sigma points around x1 
% X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
[z1,Z1,P2,Z2]=ut_knwn_inpt_nam_exp(hmeas,X,Wm,Wc,Lx,Lq);  %unscented transformation of measurments
if meas == 0 % if range (and bearing) are NaN skip correction step
    x=x1;                              %state update
    P=P1;                                %covariance update    
else
    % minimize phase difference (z(4))
    ph_err = z(4)-z1(4);
    while ph_err > pi
        z(4) = z(4) - 2*pi;
        ph_err = z(4)-z1(4);
    end
    while ph_err < -pi
        z(4) = z(4) + 2*pi;
        ph_err = z(4)-z1(4);
    end
    P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance
    K=P12*inv(P2);
    x=x1+K*(z-z1);                              %state update
    P=P1-K*P12';                                %covariance update    
end