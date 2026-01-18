% function [R_RSMA_ct,EHR_ct] = RSMA_fix(K,J,N,P_dBm_max,Emin)
function [SumRate,Loc_Conv,W] = RSMA_fix_v2(K,J,N,P_dBm_max,Emin,Loc_IDR,Loc_EHR,D_x,D_y,d,sigma2)

%解广播情况下的多波导多用户最大化最小速率问题
%天线位置固定到最近用户的位置
%% parameter
% P_dBm_max = 0;
P_max = 10.^((P_dBm_max-30)./10);

%% 
%初始化天线位置
lambda=3e8/28e9;
Loc_Conv=[zeros(1,N); D_y/2+((0:N-1)+0.5)*lambda/2; d * ones(1,N)];

SumRate = -inf;
W = NaN;

%天线到用户的信道：h_IDR[N,K],h_EHR[N,J]
h_IDR = channel_Conv(Loc_Conv, Loc_IDR);
h_EHR = channel_Conv(Loc_Conv, Loc_EHR);
[W_p, wc,SumRate] = rsma_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin);
W = [wc,W_p];
end