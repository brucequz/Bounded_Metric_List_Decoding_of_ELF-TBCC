clc
set(0,'DefaultTextFontName','Times','DefaultTextFontSize',14,...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
    'DefaultLineLineWidth',2,'DefaultLineMarkerSize',7.75);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

marksize = 1:10;

Colors = [...  
[0 0.4470 0.7410]           % blue
[0.8500 0.3250 0.0980]      % orange
[0.9290 0.6940 0.1250]      % yellow  
[0.4940 0.1840 0.5560]  	% purple      
[0.4660, 0.6740, 0.1880]	% green
[0.3010, 0.7450, 0.9330]	% light blue 
[0.6350 0.0780 0.1840]      % burgundy
[0.75, 0, 0.75]             % magenta
[0.5, 0, 0]                 % dark red
[0, 0.5, 0]                 % dark green
[0.5, 0.5, 0]               % barf green
[0.25, 0.25, 0.25]          % gray
[0.8, 0, 0]];              % red


darkBlue = Colors(1,:);
orange = Colors(2,:);
yellow = Colors(3,:);
purple = Colors(4,:);
green = Colors(5,:);
lightBlue = Colors(6,:);
burgundy = Colors(7,:);
magenta = Colors(8,:);
darkRed = Colors(9,:);
darkGreen = Colors(10,:);
barfGreen = Colors(11,:);
gray = Colors(12,:);
red = Colors(13,:);

CRC64 = [
    0x3u16  0x5u16  0xFu16  0x11u16 0x33u16 0x7Fu16 0xFFu16 0x1ABu16 0x301u16 0x4F5u16 0x9AFu16 0x1565u16];
CRC76 = [
    0x3u16  0x5u16  0xFu16  0x11u16 0x33u16 0x55u16 0x81u16 0x195u16 0x325u16 0x53Du16 0xE0Du16 0x1565u16];

N = 64;
K=64;
nu = 8;
m=12;
Eb=0;
SNRdB = 1:0.1:6;
SNRlen = length(SNRdB);
SNR = power(10,SNRdB/10);
W = exp(-SNR/2);
Wmax = 3*N-1;

Aw_name = sprintf("AwCRCandTBCCpunctV%dm%dK%d", nu,m,K);
Aw_file = sprintf('DATA/%s.mat', Aw_name);

Aw1_name = sprintf("Aw1CRCandTBCCpunctV%dm%dK%d", nu,m,K);
Aw1_file = sprintf('DATA/%s.mat', Aw1_name);

DataPunct = load(Aw1_file);
AwPunct = DataPunct.Aw1';
AwPunct(1:9)= AwPunct(10);

Data14 = load('DATA/AwV14');
Y = Data14.Aw(2:end);


Wpow = 1:Wmax;
Am = power(W, Wpow');
Aw14 = Y*Am;

dmin14 = 0;
dminPunct = 15;
for i = 1:3*N-1
    if Y(i)
        dmin14 = i;
        break
    end
end

Qcoef14 = qfunc(sqrt(dmin14*SNR));

Ecoef14 = exp(dmin14*SNR/2);
Tbound14 = Qcoef14.*Ecoef14.*Aw14;


QcoefPunct = qfunc(sqrt(dminPunct*SNR));
EcoefPunct = exp(dminPunct*SNR/2);
TboundPunct = QcoefPunct.*EcoefPunct.*AwPunct;


X = [-1 1];
pX = 0.5*ones(1,2);
% K=96;
K=64;
RCU0 = RCU_warpper(K, 2*(K), pX, X, SNRdB);
K=64;





thd = 1e-6;
logthd=-6;
for r = 1:1

    
        
    for j = 1:SNRlen-1
        if TboundPunct(j+1) < thd
%             RCUthd(m) = SNRdB(j+1);
            break;
        end
    end
   
    v1 = log10(TboundPunct(j));
    v2 = log10(TboundPunct(j+1));
        
    sv = (v2-v1)/(RCU0(j+1,1+Eb)-RCU0(j,1+Eb));
    dBpunct = RCU0(j,1+Eb)+ (logthd-v1)/sv;
    
    
        
    for j = 1:SNRlen-1
        if Tbound14(j+1) < thd
%             RCUthd(m) = SNRdB(j+1);
            break;
        end
    end
   
    
    v1 = log10(Tbound14(j));
    v2 = log10(Tbound14(j+1));
    
        
    sv = (v2-v1)/(RCU0(j+1,1+Eb)-RCU0(j,1+Eb));
    dB14 = RCU0(j,1+Eb)+ (logthd-v1)/sv;
    
    
        
    for j = 1:SNRlen-1
        if RCU0(j+1, 4) < thd
%             RCUthd(m) = SNRdB(j+1);
            break;
        end
    end
   
    v1 = log10(RCU0(j, 4));
    v2 = log10(RCU0(j+1, 4));
    
    
    sv = (v2-v1)/(RCU0(j+1,1+Eb)-RCU0(j,1+Eb));
    dBRCU0 = RCU0(j,1+Eb)+ (logthd-v1)/sv;
    

    
    
    fprintf("RCU cross: %.10f\n", dBRCU0);
    fprintf("DSU v=%d No ELF EbN0 cross %.10f, gap to RCU: %.10f\n", 14, dB14,  dB14-dBRCU0);
    fprintf("DSU m: %d: 0x%X: EbN0 cross: %.10f, gap to RCU: %.10f\n", m, CRC64(m), dBpunct, dBpunct-dBRCU0);
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

EbNo = [1,1.5,2,2.5,2.75,3,3.25,3.5, 3.75];
m12ELF = [0.1116,0.0345,0.00745793,0.00126769,0.000368773,1.0125E-04,2.65e-05,5.55E-06,1.2714e-06];

EbNoAPLVD = [1,1.5,2,2.5,3,3.25, 3.5];
APLVD = [1.12E-01, 3.14E-02, 6.67E-03, 8.43E-04, 9.43E-05, 1.81E-05, 4.38e-6];

lab = sprintf('DSU $\\nu=%2d$, $m=12$, $d_{\\min}=%d$',8, dminPunct);
semilogy( RCU0(:,1+Eb), TboundPunct,'-o', 'Color', darkBlue, 'MarkerSize', marksize(5),'LineWidth', 1.9,'DisplayName', lab);hold on


% lab = sprintf('DSU Bound CRC 0x%X dmin: %d', CRCs(minDidx), dmins(minDidx));
lab = sprintf('DSU $\\nu=14$, No ELF, $d_{\\min}=%d$', dmin14);
semilogy( RCU0(:,1+Eb), Tbound14,'-', 'Color', red, 'MarkerSize', marksize(8),'LineWidth', 1.9,'DisplayName', lab);hold on



lab = sprintf('RCU $(%d, %d)$', 2*(K), K);
semilogy( RCU0(:,1+Eb), RCU0(:,4),'--','Color', 'k', 'MarkerSize', marksize(4),'LineWidth',  1.9,'DisplayName', lab);hold on

lab = sprintf('SIM $\\nu=%2d$, $m=12$, $d_{\\min}=%d$',8, dminPunct);
semilogy( EbNo, m12ELF,'-D','Color', darkBlue, 'MarkerSize', marksize(5),'LineWidth',  1.9,'DisplayName', lab);hold on



lab = sprintf('SIM $\\nu=14$, No ELF, $d_{\\min}=%d$', dmin14);
semilogy( EbNoAPLVD, APLVD,'-S','Color', red, 'MarkerSize', marksize(5),'LineWidth',  1.9,'DisplayName', lab);hold on


grid on;
set(gca,'ycolor','k') 
% legend('location', 'northeast');

legend('location', 'southwest');

xlim([min(SNRdB) 8]);
xlim([1 6]);
ylim([1e-10 1e-1]);

yvals = -1:-1:-10;
yticklab = sprintf('$10^{%d}$\n', yvals);
yticks(10.^fliplr(yvals));
ylim([1e-12 1e-2]);
ylim([1e-10 1e-1]);

% set(gcf,'position',[600, 600, 800,600])
set(gcf,'position',[300, 400, 600,420])
if Eb
    
    xlabel('$E_b/N_0$(dB)','interpreter','latex');
    ylim([1e-6 1e-1]);
    xlim([1 4]);
else
    xlabel('SNR (dB)','interpreter','latex');
    ylim([1e-6 1e-1]);
    xlim([1 4]);
end
ylabel('Codeword Error Rate','interpreter','latex');

