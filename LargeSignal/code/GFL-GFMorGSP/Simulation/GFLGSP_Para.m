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

%% Rated line impedance1
Lg = Xg;
Rg = Rg;
L1 = X1;
R1 = R1;
L2 = X2;
R2 = R2;
position = position;
Rf=Rf;



%% AC filter parameters
Lf = 0.2;
Cf = 0.01;
Lc = 1e-9;

%% Grid-following inverter1
% PLL1
w_pll1 = k1p;   % (rad/s)
w_tau1 = 1e3 *2*pi;  % (rad/s)
kp_pll1= w_pll1;
ki_pll1= 0;

% Current loop
w_i_GFL1 = 1000 *2*pi;    % (rad/s)


%% Grid-following inverter2 - GSP
% PLL2
w_pll2 = k2p;   % (rad/s)
w_tau2 = 1e3 *2*pi;  % (rad/s)
kp_pll2= w_pll2;
ki_pll2= 0;

% Current loop
w_i_GFL2 = 1000 *2*pi;    % (rad/s)

%Iq droop control
w_tvc = 25*2*pi;
Kq = Kq;



%% voltage sag1
%voltage sag
t0_sag=10;
dt_sag=0.01;
v_sag=0.1;

%phase jump
t0_jump=10;
dt_jump=2;
phase_jump = -1.2;

%line cutting
t0 = 2.5;
tc1 = t_c-0.01;
tc2 = t_c;

