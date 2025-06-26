clc

%% Fundamental parameters            
Fs = 1e5;
Ts = 1/Fs;

%% Base values
Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

%% AC filter parameters
Lf = 0.05;
Cf = 0.01; %ideal for no Cf 
Lc = 1e-9;

% Q compensation
delta0 = asin(2*Xg*Pref/Ug^2)/2; 
V_inv = Ug*cos(delta0);
Iq0 = V_inv*Cf+Iq00;
Id0 = Ug*sin(delta0)/Xg;

%% Rated line impedance1
switch fault_type
    case "line_cut"
        Lgg = Xgg;
        Rgg = Rgg;
        Lg1 = X1;
        Rg1 = R1;
    otherwise
        Rg11 = Rg;
        Lg11 = Xg;
end


%% Grid-following inverter1
% PLL1
w_tau = 1e3 *2*pi;  % (rad/s)
kp_pll;
ki_pll;

% Current loop
w_i_GFL = 1200 *2*pi;    % (rad/s) 1000*2*pi

% DC-link
Y_dc;
Vdc_ref;
kp_v_dc;
ki_v_dc;

% ac-link
w_tvc = w_tvc;
m_v = m_v;
Vac_ref = Vac_ref ;

%% fault
% voltage sag1
switch fault_type
    case "voltage_sag"
    %voltage sag
    t_sim_start = 10;
    t0_sag = t_sim_start +t_start;
    dt_sag = t_c;
    v_sag= Ug_fault;
    case "line_cut"
    %line cutting
    case "power_jump"
     t_sim_start = 6;
     t0_jump = t_sim_start +t_start;
    case "phase_jump"
     t_sim_start = 6;
     t0_jump = t_sim_start +t_start;
end




