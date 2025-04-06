% This function generates the paramters at workspace for simulink model.

% Author(s): Yifan Zhang

%%
%clear all
%close all
clc

%% Fundamental parameters
Fs = 100;             % (kHz)
Fs = Fs*1e3;
Ts = 1/Fs;

%% Base values
Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;



%% Rated line impedance1
Lg1 = X1;
Rg1 = 0.1e-3;
Lg2 = X2;
Rg2 = 0.1e-3;
Lg0 = Xg;
Rg0 = 0.1e-3;
position = position;
Rf=Rf;



%% Grid-following inverter1
% PLL1
w_pll1 = 10 *2*pi;   % (rad/s)
w_tau1 = 1e3 *2*pi;  % (rad/s)
kp_pll1= w_pll1;
ki_pll1= 0;

% Current loop
w_i_GFL1 = 4000 *2*pi;    % (rad/s)



% AC filter parameters
Lf = 0.2;
Cf = 0.001;
Lc = 1e-9;

%% Grid-following inverter2
% PLL2
w_pll2 = 10 *2*pi;   % (rad/s)
w_tau2 = 1e3 *2*pi;  % (rad/s)
kp_pll2= w_pll2;
ki_pll2= 0;

% Current loop
w_i_GFL2 = 4000 *2*pi;    % (rad/s)

%% Fault

%voltage sag
t0_sag=10;
dt_sag=0.01;
v_sag=0.1;

%phase jump
t0_jump=10;
dt_jump=2;
phase_jump = -1.2;

%line cutting
t0 = 2;
tc1 = t_c;
tc2 = t_c;


