clear
clc
close all

figure
x = [1 2 3 4 5 6 7 8 9 10];
vals = [19 17 13 11 6 3 2 1 1 1];
b = bar(x,vals, 'r');
set(gca,'xticklabel',{'0','50','62', '75', '87', '100', '125', '150', '175', '200'})

xlabel('Transorbital Carrier Frequencies (Hz)','FontSize',12)
ylabel('Number of Subjects Who Detected Phosphenes','FontSize',12)