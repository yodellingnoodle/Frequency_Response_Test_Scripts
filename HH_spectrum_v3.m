% Hodgkin Huxley model 
% Numerically integrated with an Exponential Euler Scheme 

%Units:
%voltage is in millivolts (mV)
%current is in microamperes (uA)
%capacitance is in microfarads (uF)
%conductance is in millisiemens (mS)
%area is in centimeters squared (cm^2)
%time is in milliseconds (ms)
clear all
close all

Transfer_fnc=1 ; %% create transfer function plot
f0=10;%% create transfer function plot -initial frequency to scan
freq_step=20; %% frequency step
number_of_freq=50;%% number of  frequency to scan on
threshold=0; %% Threshold of the spike
% Area of cell
A = 1; %cm^2
%membrane capacitance per unit area:
C = 1.0;      % (uF/cm^2)
%max possible Na+ conductance per unit area:
gNabar = 120; % (mS/cm^2)
%max possible K+ conductance per unit area:
gKbar = 36;   % (mS/cm^2)
%leakage conductance per unit area:
gLbar = 0.3;  % (mS/cm^2)
%Na+ equilibrium potential:
ENa = 45;   %   (mV)
%K+ equilibrium potential:
EK = -82;   %   (mV)
%leakage channel reversal potential:
EL = -59;   %   (mV)

%initialize time step and experiment duration:
dt = 0.1;     % time step duration (ms)
tmax = 1000;    %duration of experiment (ms)

%total number of time steps in the experiment:
niter = ceil(tmax/dt);

%initialize arrays that hold data for plotting:
m_plot = zeros(1,niter);
h_plot = zeros(1,niter);
n_plot = zeros(1,niter);
V_plot = zeros(1,niter);
C_plot = zeros(1,niter);
t_plot = (0:niter-1)*dt;% time vector in ms
vstart = -70;


% Current injected
Ie = zeros(1, niter);

% square pulse of current starts at 30 ms and ends at 40 ms
%Ie(t_plot>= 30 & t_plot<= 40) = 10; % (uA/cm^2)

% square burst of  pulses of currents
%Ie((t_plot>= 30 & t_plot<= 31)|(t_plot>= 32 & t_plot<= 33)|(t_plot>= 34 & t_plot<= 35)) = 2; % (uA/cm^2)

% Single-tone stimulus -continous time
%Ie_single_tones =2*(sin(2*pi*50.*t_plot*10^-3));% Create two tones
%Ie=Ie_single_tones;

% Two-tone stimulus -continous time
%Ie_two_tones =20*(sin(2*pi*1050.*t_plot*10^-3)+sin(2*pi*1000.*t_plot*10^-3));% Create two tones
%Ie=Ie_two_tones;

% Multi-tone stimulus -continous time
%Ie_multi_tones =10*(sin(2*pi*900.*t_plot*10^-3)+sin(2*pi*950.*t_plot*10^-3)+sin(2*pi*800.*t_plot*10^-3)+sin(2*pi*850.*t_plot*10^-3)+sin(2*pi*1000.*t_plot*10^-3)+sin(2*pi*1050.*t_plot*10^-3)+sin(2*pi*500.*t_plot*10^-3)+sin(2*pi*550.*t_plot*10^-3));% Create two tones
%Ie=Ie_multi_tones;

% Two-tone stimulus -Burst two tones
%Ie_two_tones(t_plot< 30 | t_plot>= 40) = 0;
%Ie=Ie_two_tones;

% Two-tone stimulus superimposed on square stimulus
%Ie=(Ie_two_tones+Ie)*0.5;
if  Transfer_fnc==1
for i=1:number_of_freq
   
    f=f0+(i-1)*freq_step;
    Ie_single_tones =10*(sin(2*pi*f*t_plot*10^-3));% Create single tones
    Ie=Ie_single_tones;
    %voltage just at t=0: initial condition
    V_plot(1) = vstart;
    % In fact we are assuming V was at -70 before experiment to set m,h,n to 
    % their steady state values at -70

    % minf, hinf, ninf at -70  for initial condition
    m_plot(1) = alpham(V_plot(1))/(alpham(V_plot(1))+betam(V_plot(1)));
    h_plot(1) = alphah(V_plot(1))/(alphah(V_plot(1))+betah(V_plot(1)));
    n_plot(1) = alphan(V_plot(1))/(alphan(V_plot(1))+betan(V_plot(1)));

% Main for loop to numerically integrate the HH model 
 for k = 1: niter-1
  % taus for m,h,n
    tau_m = 1 /(alpham(V_plot(k))+betam(V_plot(k)));
    tau_h = 1 /(alphah(V_plot(k))+betah(V_plot(k)));
    tau_n = 1 /(alphan(V_plot(k))+betan(V_plot(k)));
  
  % Steady state values for m, h, n
    m_inf = alpham(V_plot(k))/(alpham(V_plot(k))+betam(V_plot(k)));
    h_inf = alphah(V_plot(k))/(alphah(V_plot(k))+betah(V_plot(k)));
    n_inf = alphan(V_plot(k))/(alphan(V_plot(k))+betan(V_plot(k)));  
  
  % Update m, h, and n using the Exponential-Euler method
    m_plot(k+1) = m_inf+(m_plot(k)-m_inf)*exp(-dt/tau_m);
    h_plot(k+1) = h_inf+(h_plot(k)-h_inf)*exp(-dt/tau_h);
    n_plot(k+1) = n_inf+(n_plot(k)-n_inf)*exp(-dt/tau_n);
  
  % Update conductances
    gNa = gNabar*(m_plot(k+1)^3)*h_plot(k+1);    %sodium conductance
    gK = gKbar*(n_plot(k+1)^4);                  %potassium conductance
    g = gNa+gK+gLbar;                            %total conductance
    gE = gNa*ENa+gK*EK+gLbar*EL;                 %gE=g*E
    
  % update vinf
    V_inf = (gE + Ie(k)/A) / g;
  % update tauv
%     
%     Cnl = C*(1+0.0002*(V_plot(k)-vstart).^2);
%     C_plot(k)=Cnl;
%     %C_plot(k)=C+dC*(1-exp(-dt/tau_V));
    tau_V = C/g;
  % exponential euler for updating membrane potential
    V_plot(k+1) = V_inf + (V_plot(k)-V_inf)*exp(-dt/tau_V);
end
%%%%%%%%%%%%%%[pkt,lct] = findpeaks(V_plot,t_plot,'MinPeakHeight',0);
%%%%%%%%%%%%%%%%count(i)=numel(pkt);
end
end

%%% Ie , Menbrane Voltage, Gating  plot
figure (1)
subplot(3,1,1)
plot(t_plot, V_plot, 'linewidth',2); ylim([-110 40]);legend('Vm');
%xlabel('time(ms)', 'fontsize', 20)
ylabel('V_m (mV)', 'fontsize', 20)
subplot(3,1,2)
plot(t_plot,Ie, 'LineWidth',2)
%xlabel('time(ms)', 'fontsize', 20)
ylabel('I_e (\muA/cm^2)', 'fontsize', 20)
gating = [m_plot;h_plot;n_plot];
subplot(3,1,3)
plot(t_plot, gating, 'linewidth',2);legend('m', 'h', 'n');
xlabel('time(msec)', 'fontsize', 20)

%%% Conductances plot
figure(2)
p1 = plot(t_plot,gKbar*n_plot.^4,'m','linewidth',2);
hold on
p2 = plot(t_plot,gNabar*(m_plot.^3).*h_plot,'g','linewidth',2);
legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
ylabel('Conductance')
xlabel('time (ms)')
title('Conductance for Potassium and Sodium Ions')

%%% Spectrum plot - Output / input
figure (3)
psd = abs(fft(Ie)).^2;
L=length(psd);
Fs=1/(dt*10^-3);
f=Fs/L *(1:L/2+1);
W1=psd(1:L/2+1);
W1(2:end-1) = 2*W1(2:end-1);
%plot(10*log10(W1))
%hold on
psd = abs(fft(V_plot)).^2;
W2=psd(1:L/2+1);
W2(2:end-1) = 2*W2(2:end-1);
HH_Spectrum_plot =[10*log10(W1);10*log10(W2)]

plot(f,HH_Spectrum_plot);legend('Ie input ', 'Menbrane potential ')
title('Single-Sided Spectrum of HH neuron to stimulation')
xlabel('f (Hz)')
ylabel('Ie & Action potential [dB]')


%%% Spectrum plot - gating variables

figure(4)
psd_nplot = abs(fft(n_plot)).^2;
W1=psd_nplot(1:L/2+1);
W1(2:end-1) = 2*W1(2:end-1);
psd_mplot = abs(fft(m_plot)).^2;
W2=psd_mplot(1:L/2+1);
W2(2:end-1) = 2*W2(2:end-1);
Gating_variable_Spectrum = [10*log10(W1);10*log10(W2)];

plot(f,Gating_variable_Spectrum);
legend('Potassium gating variable spectrum', 'Sodium gating variable spectrum')
title('Single-Sided Spectrum of gating variable')
xlabel('f (Hz)')
ylabel('Gating variable spectrum [dB]')


figure(5)
plot(t_plot, C_plot, 'linewidth',2);
xlabel('time(msec)', 'fontsize', 20)
legend('Capacitance over time')
title('Non-Linear Capacitance voltage dependend')
xlabel('time(msec)')
ylabel('Capacitance (uF)')

figure(6)
plot(V_plot, C_plot, 'linewidth',2);
xlabel('Voltage(mV)', 'fontsize', 20)
legend('Capacitance over time')
title('Non-Linear Capacitance voltage dependend')
xlabel('Voltage(mV)')
ylabel('Capacitance (uF)')

figure(7)
plot(V_plot, tau_V, 'linewidth',2);
xlabel('Voltage(mV)', 'fontsize', 20)
legend('Tau over time')
title('Tau voltage dependend')
xlabel('Voltage(mV)')
ylabel('Tau (msec)')

%%% Spiking rate transfer function  plot
if  Transfer_fnc==1
figure(8)
i=1:number_of_freq;
p1 = plot((f0+(i-1)*freq_step),count,'m','linewidth',2);
legend([p1], 'Spiking rate as function of stimulation frequency')
ylabel('rate')
xlabel('freq (Hz)')
title('Spiking rate transfer function  plot')
end
