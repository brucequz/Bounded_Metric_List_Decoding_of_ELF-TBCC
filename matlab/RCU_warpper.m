function FER = RCU_warpper(k, n, pX, X, snr)
%--------------------------------------------------------------------------
%       FER = RCU_warpper(k, n, pX, X, snr)
% This function is a wrapper function of  Fabian's RCU bound
% Input     k       --      input bits 
%           n       --      block length (real channel)
%           pX      --      channel input distirbuion
%           X       --      channel input set
%           snr     --      SNR dB:
% Output    FER     --      a matrix that contains SNR/EbN0/EsN0 and FER
%                           column-1 SNR
%                           column-2 EbN0
%                           column-3 EsN0
%                           column-4 FER
%
%                                                   Linfang Wang
%--------------------------------------------------------------------------

rate = k/n;
FER =zeros(length(snr),4);
FER(:,1)=snr.';
FER(:,2)=FER(:,1)-10*log10(2*rate); %EB/N0
FER(:,3)=FER(:,2)+10*log10(rate); %ES/N0
%le input channel inputs e
ave_P = 0;
for ii = 1 : length(X)
    ave_P=ave_P+pX(ii)*X(ii)^2;
end
X = X/sqrt(ave_P);
for ii = 1:length(snr)
    sigma2=10^(-snr(ii)/10);
    FER(ii,4)=rcu(n,rate, pX.' , X.' , sigma2); 
end
end