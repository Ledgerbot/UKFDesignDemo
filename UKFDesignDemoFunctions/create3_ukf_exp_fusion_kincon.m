% This file simulates the UKF on the create 3 robot. 
% The outputs are the state estimate (xh_k) and the associated covariance
% (Px_k)
% The inputs are 
% xh0    = initial state estimate
% Px0    = initial state covariance
% Rv     = process noise covariance
% Rn     = measurement noise covariance
% yMeas  = measured data (includes multiple sensor outputs by default)
% uest   = 1 to estimate speed, 0 if speed is known
% n_sens = # of independent sensors
% p_msng = percentage of missing measurements

function agents = create3_ukf_exp_fusion_kincon(xh0,Px0,Rv,Rn,Xa,Ya,Yaw,YawRt,U,Y,Ts,Na,umax,broadcast)
% Uses code developed by Yi Cao 
% https://www.mathworks.com/matlabcentral/fileexchange/18217-learning-the-unscented-kalman-filter
% Modified by R. O'Brien to incorporate non-additive process noise (5 July 2023)
% Modified by R. O'Brien to incorporate consensus (17 June 2024)
% Modified by R. O'Brien to extend to N agents (24 June 2024)
% Modified by R. O'Brien for use with experimental data (27 July 2024)

%% UKF parameters

% # of measurements
m = 4;

% # of states
n = size(xh0,1);

% process noise input size
q = n;

% Total dimension for UKF
N = n + q + m;

% Assume 100% broadcast

%% Initialize variables
% Each agent has estimate/covariance for self and other agents as well as
% consensus target estimate/covariance

Nsteps        = size(Xa,1);

% Initialize agents 
for i = 1:Na % agents    
    for j = 1:Na % virtual agents for full federated process      
        agents(i,j).xh_k        = zeros(n,Nsteps);
        agents(i,j).xh_k(:,1)   = xh0(:,j);
        agents(i,j).Px_k        = zeros(n,n,Nsteps);
        agents(i,j).Px_k(:,:,1) = Px0;
        agents(i,j).xhm         = zeros(4,1);
        agents(i,j).Pxm         = zeros(4,4);
    end
end


%% federated fusion
aset = 1:Na;
ro = 2; % bearing only
for k = 1:Nsteps-1             
    % Update the self-estimate of the pose and victim state for each agent
    
    for i = 1:Na % loop thru agents
        % previous (sample k) self-estimate and self-covariance
        xhp = agents(i,i).xh_k(:,k);
        Pxp = squeeze(agents(i,i).Px_k(:,:,k));
        Uk = [YawRt(k,i);U(k,i)]; % known input for agent i
        if ~isnan(Y(k,i))
            meas = 1;
            yk = [Xa(k,i);Ya(k,i);Yaw(k,i);Y(k,i)]; % measurement is agent's state and bearing to target        
            [x,P,xm,Pm] = ukf_exp(@create3_exp,xhp,Pxp,@c3_meas_exp,yk,Rv,Rn,Uk,meas,ro,Ts);            
        else
            meas = 0;
            yk = zeros(4,1); % measurement is agent's state and bearing to target            
            [x,P,xm,Pm] = ukf_exp(@create3_exp,xhp,Pxp,@c3_meas_exp,yk,Rv,Rn,Uk,meas,ro,Ts);
        end
        if x(7)< 0
            x(7) = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % limit change in estimated target position using kinematic
        % constraint on forward velocity
        ut_est = norm(x(4:5)-xhp(4:5))/Ts;
        if ut_est > umax
            x(4:5) = xhp(4:5) + (umax/ut_est)*(x(4:5)-xhp(4:5));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store current estimate/covariance
        agents(i,i).xh_k(:,k+1)   = x;
        agents(i,i).Px_k(:,:,k+1) = P;
        agents(i,i).xhm           = xm(4:7);
        agents(i,i).Pxm           = Pm(4:7,4:7);
    end
    
    % Use broadcast information to update estimate of other agents'
    % estimates
    for i1 = 1:Na % all agents
        jset = setdiff(aset,i1);
        for j = 1:Na-1 % update estimates of other agents but not self (updated above)
            j1 = jset(j);           
            if broadcast(j1,k) == 1  % agent j broadcasts its self state/covariance estimate
                x = agents(j1,j1).xh_k(:,k+1);
                P = agents(j1,j1).Px_k(:,:,k+1);
            else % agent j did not broadcast
                % use agent i's previous (sample k) estimate and covariance for agent j
                xhp = agents(i,j1).xh_k(:,k);
                Pxp = squeeze(agents(i,j1).Px_k(:,:,k));
                Uk = zeros(2,1); % (agent j's input is unknown) squeeze(AI(k,:,j1))'; % known input for agent j
                meas = 0;            
                yk = zeros(4,1); % no measurement received
                [x,P,xm,Pm] = ukf_exp(@create3_exp,xhp,Pxp,@c3_meas_exp,yk,Rv,Rn,Uk,meas,ro,Ts);
                if x(7)< 0
                    x(7) = 0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % limit change in estimated target position using kinematic
                % constraint on forward velocity
                ut_est = norm(x(4:5)-xhp(4:5))/Ts;
                if ut_est > umax
                    x(4:5) = xhp(4:5) + (umax/ut_est)*(x(4:5)-xhp(4:5));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end % meas loop
            % Store current estimate/covariance
            agents(i1,j1).xh_k(:,k+1)   = x;
            agents(i1,j1).Px_k(:,:,k+1) = P;
            agents(i1,j1).xhm           = xm(4:7);
            agents(i1,j1).Pxm           = Pm(4:7,4:7);
        end % j loop
    end
    % Fusion: Use conformity if only one agent receives a measurement
    for i2 = 1:Na
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Info fusing: federated for agent i (full)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % federated estimate/covariance across agent i's UKF for self and other agents
        % previous federated estimate/covariance (use full.agent 1 b/c all
        % agents reset to same values
        Pxp = squeeze(agents(i2,1).Px_k(:,:,k)); % all agents have same covariance after reset
        xhp = agents(i2,1).xh_k(:,k); % all agents have same estimate after reset
        meas = 0;
        yk = zeros(5,1); % no measurement received
        Uk = [YawRt(k,i);U(k,i2)]; % known input for agent i
        [x,P] = ukf_exp(@create3_exp,xhp,Pxp,@c3_meas_exp,yk,Rv,Rn,Uk,meas,ro,Ts);        
        if x(7)< 0
            x(7) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % limit change in estimated target position using kinematic
        % constraint on forward velocity
        ut_est = norm(x(4:5)-xhp(4:5))/Ts;
        if ut_est > umax
            x(4:5) = xhp(4:5) + (umax/ut_est)*(x(4:5)-xhp(4:5));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xhm = x(4:7);  
        Pxm = P(4:7,4:7);
        O_sum = zeros(4,4);
        phi_sum = zeros(4,1);
        for j2 = 1:Na
            Oj    = inv(squeeze(agents(i2,j2).Px_k(4:7,4:7,k+1)));
            Ojm   = inv(squeeze(agents(i2,j2).Pxm));
            xj    = agents(i2,j2).xh_k(4:7,k+1);
            xjm   = agents(i2,j2).xhm;
            O_sum = O_sum + Oj-Ojm;
            phi_sum = phi_sum + Oj*xj-Ojm*xjm;
        end
        Of  = inv(Pxm) + O_sum;                % federated inverse covariance
        xf  = inv(Of)*(inv(Pxm)*xhm+phi_sum);  % federated (target) state estimate

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % limit change in estimated target position using kinematic
        % constraint on forward velocity
        % if kincon
        ut_est = norm(xf(1:2)-xhp(4:5))/Ts;
        if ut_est > umax
            xf(1:2) = xhp(4:5) + (umax/ut_est)*(xf(1:2)-xhp(4:5));
        end
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        % reset all virtual agents to the federated target estimate and
        % covariance values
        for j3 = 1:Na
            agents(i2,j3).xh_k(4:7,k+1)      = xf;
            agents(i2,j3).Px_k(4:7,4:7,k+1)  = inv(Of);            
        end

        % Outlier ID and rejection code (TBD)

    end % agents loop
end % k loop