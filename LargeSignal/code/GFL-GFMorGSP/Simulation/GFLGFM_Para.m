% This function generates the paramters at workspace for simulink model.

% Author(s): Yifan Zhang

%%
%clear all
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

%% AC filter parameters
Lf1 = 0.1;
Cf1 = 0.002;
Lc1 = 1e-9;

Lf2 = 0.1;
Cf2 = 0.2;
Lc2 = 1e-9;

%% Grid-following inverter
Id = I1d;
% PLL
w_pll1 = k1p;     % (rad/s)
w_tau1 = 1e3*2*pi;    % (rad/s)
kp_pll1= w_pll1;
ki_pll1= 0;

% Current loop
w_i_GFL1 = 2500 *2*pi;    % (rad/s)

%% Grid-forming inverter
% Droop
Pm = Pref;
M = Kgfm/Wbase;        % (pu)
w_droop = 15*2*pi;   %10 1 0

% Current loop
w_i_GFM = 2500*2*pi;

% Voltage loop
w_v_GFM = 200 *2*pi;
Scale_ki_v = 60;


%% Grid parameters
Lg = Xg;
Rg = Rg;
L1 = X1;
R1 = R1;
L2 = X2;
R2 = R2;


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
tc2 = t_c2;
