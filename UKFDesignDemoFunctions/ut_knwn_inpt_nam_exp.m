function [y,Y,P,Y1]=ut_knwn_inpt_nam_exp(f,X,Wm,Wc,n,nq)
%Unscented Transformation (non additive noise)
%Input:
%        f: nonlinear map (state, process noise)
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%       nq: numer of process noise inputs 
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations
L1=size(X,2);
LX= size(X,1);

nr = LX-n-nq-1;

y=zeros(nr,1);
Y=zeros(nr,L1);
for k=1:L1                   
    Y(:,k)=f(X(1:n,k),X(n+1+nq:LX,k));       
    y=y+Wm(k)*Y(:,k);       
end
Y1=Y-y(:,ones(1,L1));
P=Y1*diag(Wc)*Y1'; % + 0.01*eye(nr);          