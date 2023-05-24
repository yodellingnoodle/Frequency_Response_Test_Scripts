clear
clc
close all

figure
x = [1 2 3 4 5 6 7 8 9 10];
vals = [18 17 17 16 12 10 6 2 2 5];
b = bar(x,vals);
set(gca,'xticklabel',{'0','50','75', '100', '250', '500', '750', '1000', '2500', '5000'})
ylim([0, 20])

xlabel('Forearm Carrier Frequencies (Hz)','FontSize',12)
ylabel('Number of Subjects Who Detected Sensation','FontSize',12)