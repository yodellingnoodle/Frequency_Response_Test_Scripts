close all
clc
clear

figure
hold on
x = [0 50 75 100 250 500 750 1000 2500 5000];
y = [0.631 0.933 1.035 1.002 1.394 1.560 1.732 1.832 1.921 1.848]; 
err = [0.116 0.132 0.157 0.137 0.131 0.133 0.117 0.118 0.059 0.123];
plot(x,y,"bo","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
[p1,errorest1] = polyfit(x,y,4);
[y2,delta1] = polyval(p1,x,errorest1);
plot(x,y2,"b-")
plot(x,y2+(2*delta1),'b:')
plot(x,y2-(2*delta1),'b:')
arm_correlation = corrcoef(x,y);

xlabel('Forearm Carrier Frequencies (Hz)','FontSize',12)
ylabel('Detection Threshold (mA)','FontSize',12)
xlim([0, 5300])
ylim([0, 2.3])

