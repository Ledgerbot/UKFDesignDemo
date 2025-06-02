function [y,Y,P,Y1]=ut_knwn_inpt_exp(f,X,Wm,Wc,n,S,Uk,Ts)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations
N1 = size(X,2);
y  = zeros(n,1);
Y  = zeros(n,N1);
nu = length(Uk);
for k=1:N1                   
    Y(:,k) = f(X(1:n,k),Uk,zeros(nu,1),Ts);       % process noise inputs set to zero b/c additive term (S)
    y      = y+Wm(k)*Y(:,k);       
end
Y1 = Y-y(:,ones(1,N1));
P  = Y1*diag(Wc)*Y1'+ S;          