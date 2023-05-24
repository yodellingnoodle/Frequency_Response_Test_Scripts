close all
clc
clear

figure
hold on

a = [0 50 62 75 87 100 125 150 175 200];
b = [0.212 0.851 1.389 1.635 1.938 1.933 1.953 1.961 1.964 1.972];
rre = [0.073 0.106 0.112 0.103 0.087 0.057 0.177 0.168 0.158 0.123];
plot(a,b,"ro","MarkerSize",5,"MarkerEdgeColor","red","MarkerFaceColor",[0.65 0.85 0.90])
[p2,errorest2] = polyfit(a,b,6);
[b2,delta2] = polyval(p2,a,errorest2);
plot(a,b2,"r-")
plot(a,b2+(2*delta2),'r:')
plot(a,b2-(2*delta2),'r:')
ret_correlation = corrcoef(a,b);

xlabel('Transorbital Carrier Frequencies (Hz)','FontSize',12)
ylabel('Detection Threshold (mA)','FontSize',12)
xlim([0, 210])
ylim([0, 2.3])


