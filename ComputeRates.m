
function [sumrate,DLSumRate,ULSumRate] = ComputeRates(H_dl,H_ul,Hli,g_ul_dl,Q_dl,q_ul,NoisePower_dl,NoisePower_ul)
                        

nUser_dl = size(H_dl,2);      % The number of DL users
nUser_ul = size(H_ul,2);      % The number of UL users

nRx = size(Hli,1);           % number of receive antennas at BS on UL

%Calculate DL SumRate:
%---------------------
DLSumRate = 0;  

for iUser_dl=1:nUser_dl
    Interf_ul_dl = sum(q_ul.*abs(g_ul_dl(:,iUser_dl)).^2);
    Interf_dl = NoisePower_dl + real(H_dl(:,iUser_dl)'*Q_dl(:,:,1:nUser_dl~=iUser_dl)*H_dl(:,iUser_dl));
    SINR_dl = real(H_dl(:,iUser_dl)'*Q_dl(:,:,iUser_dl)*H_dl(:,iUser_dl))/(Interf_dl+Interf_ul_dl);
    
    %/SumRate in the DL channel:
    DLSumRate  = DLSumRate+ double(log2(1+SINR_dl));
   
 end
self_interference = Hli*sum(Q_dl,3)*Hli';
%Calculate UL SumRate: 
%---------------------
ULSumRate = 0;
for jUser_ul=1:nUser_ul

    Interf_ul = self_interference+NoisePower_ul*eye(nRx,nRx);
    Interf_ul = Interf_ul + H_ul(:,jUser_ul+1:nUser_ul)*diag(q_ul(jUser_ul+1:nUser_ul))*H_ul(:,jUser_ul+1:nUser_ul)';

    % compute SINR for each user
    SINR_ul=q_ul(jUser_ul)*real(H_ul(:,jUser_ul)'*Interf_ul^(-1)*H_ul(:,jUser_ul));
    % compute SR in the UL channel
    ULSumRate  = ULSumRate+double(log2(1+SINR_ul));
end

%Calculate Total SumRate
sumrate = DLSumRate+ULSumRate;
            

