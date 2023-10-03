clc;
clear;
close all;

data=load('mydata');

Symbols=data.Symbols;
R=data.R;

alpha=0.95;

port=PortfolioCVaR();
port=port.setScenarios(R);
port=port.setDefaultConstraints();
port=port.setProbabilityLevel(alpha);

W=port.estimateFrontier(100);
WReturn=port.estimatePortReturn(W);
WRisk=port.estimatePortRisk(W);

WD=eye(port.NumAssets);
MU=port.estimatePortReturn(WD);
SIGMA=port.estimatePortRisk(WD);

figure;
plot(WRisk,WReturn,'LineWidth',2);
hold on;
for i = 1:numel(Symbols)
    s = Symbols{i};
    mu = MU(i);
    sigma = SIGMA(i);
    plot(sigma,mu,'ro','MarkerFaceColor','r');
    text(sigma+0.001,mu,s);
end
grid on;
xlabel('Risk (CVaR)');
ylabel('Return (Mean)');
