rng(1)
%*********************%
%  SYSTEM PARAMETERS  %
%*********************%
BW = 10*1e6;      % Bandwidth of system (Hz)
No = -174;         % noise power density, dBm/Hz
NF_dl = 9;         % receiver noise figure in DL in dB
NF_ul = 5;         % receiver noise figure in UL in dB     
NoisePower_dl_dBm = No+10*log10(BW)+NF_dl;
NoisePower_dl = 10^((NoisePower_dl_dBm-30)/10);      % W/Hz noise power (DL)
NoisePower_ul_dBm = No+10*log10(BW)+NF_ul;
NoisePower_ul = 10^((NoisePower_ul_dBm-30)/10);      % W/Hz noise power (UL)

%/Self interference cancellation level:
deta_dB = -90;
deta =   10.^(deta_dB/10);


nTx = 4;          % number of transmit antennas at BS on DL
nRx = 4;           % number of receive antennas at BS on UL
% 

BS_coord = [0 0];     % coordinate of the BS (at the origin)

%/Distance b/w BS and Users:
R_c = 100;        % Radius of the users dropping are (im meters)
R_h = 10;          % Distance b/w UE and BS < R_h --> that UE don't impact on Pathloss (m)


%% positions of downlink users (specifics for 2 users)
angle_dl = [0;200];             % angle (in degree)
nUsers_dl = length(angle_dl);      % The number of DL users

radius_dl =[R_c;0.5*R_c];       % distance to the origin
% conver to rectangular coordinate
UE_dl_rec_coord = [radius_dl.*cosd(angle_dl) radius_dl.*sind(angle_dl) ];

%% coordinate of uplink users (specifics for 2 users)
angle_ul = [90;300];            % angle (in degree)
radius_ul =[0.85*R_c;0.75*R_c]; % distance to the origin
UE_ul_rec_coord = [radius_ul.*cosd(angle_ul) radius_ul.*sind(angle_ul)];
nUsers_ul = length(angle_ul);      % The number of UL users

%% Plot the positions of users
%{
plot(R_c*cos(linspace(0,2*pi,300)),R_c*sin(linspace(0,2*pi,300)),'--')
hold on
plot(UE_dl_rec_coord(:,1),UE_dl_rec_coord(:,2),'Marker','square','Color','blue',...
    'LineStyle','none','MarkerFaceColor','blue')
plot(UE_ul_rec_coord(:,1),UE_ul_rec_coord(:,2),'Marker','^','Color','k',...
    'LineStyle','none','MarkerFaceColor','none')
plot(BS_coord(1),BS_coord(2),'Marker','o','Color','red',...
    'LineStyle','none','MarkerFaceColor','red')
text(BS_coord(1),BS_coord(2)-7,'BS')
text(UE_dl_rec_coord(1,1)-15,UE_dl_rec_coord(1,2),'$D_1$','Interpreter','latex')
text(UE_dl_rec_coord(2,1)-15,UE_dl_rec_coord(2,2),'$D_2$','Interpreter','latex')
text(UE_ul_rec_coord(1,1),UE_ul_rec_coord(1,2)-10,'$U_1$','Interpreter','latex')
text(UE_ul_rec_coord(2,1),UE_ul_rec_coord(2,2)-10,'$U_2$','Interpreter','latex')
title('User Positions')
axis square;
%}

%% Power parameters
Pd_dbm = 26; % max Tx power at BS (dBm)
Pd_W = 10^((Pd_dbm-30)/10);
Pu_dbm = 23; % max Tx power at users (dBm)
Pu_W = 10^((Pu_dbm-30)/10).*ones(nUsers_ul,1);

%% Calculate distances between nodes
% distances b/w uplink user to downlink user
dist_ul_dl = zeros(nUsers_ul,nUsers_dl);
for iUser_dl=1:nUsers_ul
    for iUser_ul=1:nUsers_dl
        dist_ul_dl(iUser_dl,iUser_ul) = norm(UE_ul_rec_coord(iUser_dl,:)-UE_dl_rec_coord(iUser_ul,:));
    end
end
% distances b/w base station to downlink user
dist_BS_dlUEs = zeros(nUsers_dl,1);
for iUser_dl=1:nUsers_dl
    dist_BS_dlUEs(iUser_dl) = norm(BS_coord-UE_dl_rec_coord(iUser_dl,:));
end

% distances b/w base station to uplink user
dist_BS_ulUEs = zeros(nUsers_ul,1);
for iUser_ul=1:nUsers_ul
    dist_BS_ulUEs(iUser_ul) = norm(BS_coord-UE_ul_rec_coord(iUser_ul,:));
end

%% Path loss and shadowing
fc = 2;               % in GHz

PL_dl = zeros(nUsers_dl,1);
for iUser_dl=1:nUsers_dl
    PL_dl(iUser_dl)= 103.8+20.9*log10(dist_BS_dlUEs(iUser_dl)/1000)+20*log10(fc/2);       %Pico LOS d(km)
    PL_dl(iUser_dl)=10.^(-PL_dl(iUser_dl)/10);
end


PL_ul = zeros(nUsers_ul,1);
for iUser_ul=1:nUsers_ul
    PL_ul(iUser_ul)= 103.8+20.9*log10(dist_BS_ulUEs(iUser_ul)/1000)+20*log10(fc/2);       %Pico LOS d(km)
    PL_ul(iUser_ul)=10.^(-PL_ul(iUser_ul)/10);
end

PL_ul_dl = zeros(nUsers_ul,nUsers_dl);
for iUser_ul=1:nUsers_ul
    for iUser_dl=1:nUsers_dl
        PL_ul_dl(iUser_ul,iUser_dl)= 145.4+37.5*log10(dist_ul_dl(iUser_ul,iUser_dl)/1000);      % NLOS Pico
        PL_ul_dl(iUser_ul,iUser_dl) = 10.^(-PL_ul_dl(iUser_ul,iUser_dl)/10);
    end
end


%% Channels Generation
% DL channel
H_dl = zeros(nTx,nUsers_dl);
for n=1:nUsers_dl
    H_dl(:,n)= ((sqrt(PL_dl(n))*(randn(1,nTx)+1i*randn(1,nTx)))/sqrt(2))';
end

% UPLINK users to DOWNLINK users 
g_ul_dl = zeros(nUsers_ul,nUsers_dl);
for n=1:nUsers_ul
    for m=1:nUsers_dl
        g_ul_dl(n,m)= (sqrt(PL_ul_dl(n,m))*(randn(1,1)+1i*randn(1,1)))/sqrt(2);
    end
end
% UL channel
H_ul = zeros(nRx,nUsers_ul);
for n=1:nUsers_ul
    H_ul(:,n)= (sqrt(PL_ul(n))*(randn(nRx,1)+1i*randn(nRx,1)))/sqrt(2);
end
% Sefl Interference (Rician)
K_rician_dB = 1;
K_rician = 10.^(K_rician_dB/10);
H_ray = (randn(nRx,nTx)+1i*randn(nRx,nTx))/sqrt(2); %CN(0,1)
H_los = ones(nRx,nTx);
Hli =(sqrt(K_rician*deta/(K_rician+1))*H_los+sqrt(deta/(K_rician+1))*H_ray);

%% Data scaling
scale_dl = NoisePower_dl;
H_dl= H_dl/sqrt(scale_dl); % scale DL channels with scaling factor of scale_dl
g_ul_dl(n,m)= g_ul_dl(n,m)/sqrt(scale_dl); % scale UL-2-DL channels with scaling factor of scale_dl
EffNoisePower_dl = NoisePower_dl/scale_dl; % the effective noise power after DL channels scaling

scale_ul = NoisePower_ul; 
H_ul = H_ul/sqrt(scale_ul); % scale UL channels with scaling factor of scale_ul 
Hli = Hli/sqrt(scale_ul); % scale SI channel with scaling factor of scale_ul
EffNoisePower_ul = NoisePower_ul/scale_ul; % effective noise power after UL channels scaling
%% generate initial points
q_ul_init = Pu_W;
Q_dl_init_tmp = zeros(nTx,nTx,nUsers_dl);     % beamformer
for iUser_dl=1:nUsers_dl
    a = rand(nTx,nTx)+1i*rand(nTx,nTx);
    Q_dl_init_tmp(:,:,iUser_dl) = a*a';
end
P_tmp = real(trace(sum(Q_dl_init_tmp,3)));
Q_dl_init = Q_dl_init_tmp/(P_tmp/(Pd_W));
mypsi = zeros(nUsers_dl,1);
SINR_dl = zeros(nUsers_dl,1);
t_dl_init = zeros(nUsers_dl,1);
for iUser_dl=1:nUsers_dl
    interf_dl = EffNoisePower_dl;
    interf_dl = interf_dl + real(H_dl(:,iUser_dl)'*sum(Q_dl_init(:,:,1:nUsers_dl~=iUser_dl),3)*H_dl(:,iUser_dl));
    interf_dl = interf_dl + sum(q_ul_init.*abs(g_ul_dl(:,iUser_dl)).^2);
    SINR_dl(iUser_dl) = real((H_dl(:,iUser_dl)'*Q_dl_init(:,:,iUser_dl)*H_dl(:,iUser_dl)))/interf_dl;
    t_dl_init(iUser_dl) = 1+SINR_dl(iUser_dl);
    mypsi(iUser_dl) = (1+SINR_dl(iUser_dl))/interf_dl;
end
self_interference = EffNoisePower_ul*eye(nRx,nRx)+Hli*sum(Q_dl_init,3)*Hli';
X_init = zeros(nRx,nRx,nUsers_ul);
x_ul_init = zeros(nUsers_ul,1);
for iUser_ul=1:nUsers_ul

    %X_init(:,:,iUser_ul)= interf_ul;
    X_init(:,:,iUser_ul) = self_interference+H_ul(:,iUser_ul+1:nUsers_ul)...
        *diag(q_ul_init(iUser_ul+1:nUsers_ul))*H_ul(:,iUser_ul+1:nUsers_ul)';
    x_ul_init(iUser_ul)= sqrt(q_ul_init(iUser_ul));
end

ops=sdpsettings('solver','mosek','verbose',0); % solver options
% Define optimization variables
t_dl = sdpvar(nUsers_dl,1);
t_ul = sdpvar(nUsers_ul,1);
q_ul = sdpvar(nUsers_ul,1); % uplink power 
Q_dl = sdpvar(nTx,nTx,nUsers_dl,'hermitian','complex'); % downlink covariances
mybeta = sdpvar(nUsers_dl,1) ;
x_ul = sdpvar(nUsers_ul,1);
X = zeros(nRx,nRx,nUsers_ul);

% define (fixed) constraints:
F=[];
for iUser_dl=1:nUsers_dl
    F=[F,Q_dl(:,:,iUser_dl) >= 0];
end
F=[F,real(trace(sum(Q_dl,3))) <= Pd_W];
F=[F,q_ul >= 0];
F=[F,q_ul <= Pu_W];
F=[F,t_dl >= 1, t_ul >= 1]; % (24i)

%% Generate initial points
%% Iterative process        
obj=geomean([t_dl;t_ul]); % the objective in (24a)
MAX_ITER = 30; % max. number of iterations
lowerbound = zeros(MAX_ITER,1);
sumrate_true = zeros(MAX_ITER,1);

for iIter = 1:MAX_ITER

    F1 = [];
    for iUser_dl=1:nUsers_dl
        %/Constraint:
        RHS24b = EffNoisePower_dl; % the right-hand-side of (24b)
        RHS24b = RHS24b+real(H_dl(:,iUser_dl)'*sum(Q_dl,3)*H_dl(:,iUser_dl));
        RHS24b = RHS24b + sum(q_ul.*abs(g_ul_dl(:,iUser_dl)).^2);
        
        % Constraint (24b)
        F1=[F1,rcone([sqrt(1/(2*mypsi(iUser_dl)))*t_dl(iUser_dl);sqrt(mypsi(iUser_dl)/2)* mybeta(iUser_dl)],RHS24b,1/2)]; 

        % Constraint (24d)
        LHS24d = RHS24b - real(H_dl(:,iUser_dl)'*Q_dl(:,:,iUser_dl)*H_dl(:,iUser_dl)); % the left-hand-side of (24d)
        F1=[F1,LHS24d <= mybeta(iUser_dl)];

    end
    self_interference = EffNoisePower_ul*eye(nRx,nRx)+ Hli*sum(Q_dl,3)*Hli';

    for iUser_ul=1:nUsers_ul
        Interf_ul = self_interference +H_ul(:,iUser_ul+1:nUsers_ul)*diag(q_ul(iUser_ul+1:nUsers_ul))*H_ul(:,iUser_ul+1:nUsers_ul)';
        % Constraint (24c)
        F1=[F1,t_ul(iUser_ul)<= 1+ real(x_ul_init(iUser_ul)^2*H_ul(:,iUser_ul)'*X_init(:,:,iUser_ul)^(-1) *H_ul(:,iUser_ul))...
            + real((2*x_ul_init(iUser_ul)*H_ul(:,iUser_ul)'*X_init(:,:,iUser_ul)^(-1) *H_ul(:,iUser_ul))'*(x_ul(iUser_ul)-x_ul_init(iUser_ul)))...
            - real(trace((x_ul_init(iUser_ul)^2*(X_init(:,:,iUser_ul)^(-1))'*H_ul(:,iUser_ul)*H_ul(:,iUser_ul)'*(X_init(:,:,iUser_ul)^(-1))')'*...
            (Interf_ul-X_init(:,:,iUser_ul))))];

        %/Constraint (24e)
        F1=[F1,cone([2*x_ul(iUser_ul);1-q_ul(iUser_ul)], 1+q_ul(iUser_ul))];
    end

    diagnotics=solvesdp([F,F1],-obj,ops);
    if(diagnotics.problem==0) % successfully solved, optimal solution found
        lowerbound(iIter) = real((nUsers_ul+nUsers_dl)*log2(double(obj)));
        q_ul_init = double(q_ul);
        Q_dl_init = double(Q_dl);
        x_ul_init = double(x_ul);
        mypsi = double(t_dl./mybeta);
        self_interference = EffNoisePower_ul*eye(nRx,nRx)+ Hli*sum(Q_dl_init,3)*Hli';

        for iUser_ul=1:nUsers_ul
            X_init(:,:,iUser_ul) = self_interference+H_ul(:,iUser_ul+1:nUsers_ul)...
                *diag(q_ul_init(iUser_ul+1:nUsers_ul))*H_ul(:,iUser_ul+1:nUsers_ul)';
        end
        sumrate_true(iIter) = ComputeRates(H_dl,H_ul,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul);
    elseif ((3 <= diagnotics.problem) && (diagnotics.problem<=5))% some numerical issues, but could continue
        disp(strcat(diagnotics.info," at iteration ",num2str(iIter)))
        disp('Try using the current solution to continue')
        lowerbound(iIter) = real((nUsers_ul+nUsers_dl)*log2(double(obj)));
        q_ul_init = double(q_ul);
        Q_dl_init = double(Q_dl);
        x_ul_init = double(x_ul);
        mypsi = double(t_dl./mybeta);
        self_interference = EffNoisePower_ul*eye(nRx,nRx)+ Hli*sum(Q_dl_init,3)*Hli';

        for iUser_ul=1:nUsers_ul
            X_init(:,:,iUser_ul) = self_interference+H_ul(:,iUser_ul+1:nUsers_ul)...
                *diag(q_ul_init(iUser_ul+1:nUsers_ul))*H_ul(:,iUser_ul+1:nUsers_ul)';
        end
        sumrate_true(iIter) = ComputeRates(H_dl,H_ul,Hli,g_ul_dl,Q_dl_init,q_ul_init,EffNoisePower_dl,EffNoisePower_ul);
    else
        disp(strcat(diagnotics.info,' at iteration ',num2str(iIter)))
        sumrate_true(iIter:end) = [];
        lowerbound(iIter:end) = [];
        break
    end

end
figure
plot(lowerbound)
hold on
plot(sumrate_true)
legend('Lower bound of SumRate','True Sum Rate','Location','southeast')
saveas(gcf, '../../results/ConvergencePlot_Algorithm2.png')