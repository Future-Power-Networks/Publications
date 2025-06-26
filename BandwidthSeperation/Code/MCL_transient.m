close all;
clc;
%% parameter
%simulation
t_end =  4;

%grid-postfault
Rg = 2e-3;
Xg = 1/2.1;%2.1
Ug = 1;
Ws = 2*pi*50; 
Lg= Xg/Ws;

% initial time for boundary layer
t_0 = 1e-6;


%GFL
Id = 1;
Iq00 = 0; %negative Q (positive direction)   constant value initial value
try
    load("initial.mat");
    Iq_clear;   % for slow TVC\
    Iq_pre; % for slow TVC\
catch
    Iq_clear = Iq00;
    Iq_pre = Iq00;
end
Iq = Iq00;   % variable for function
w_pll = 16*2*pi;     % (rad/s)
kp_pll = w_pll;
ki_pll = kp_pll/4;%0;%w_pll^2/4;
w_limit = 500;

Vdc_ref = 2.5;
Y_dc = 12.5;  %12.5
C_dc = Y_dc/Ws;
w_vdc   = 2*2*pi;
kp_v_dc	= Y_dc/2*w_vdc/Ws;
zeta_dc = 1/4;
ki_v_dc	= kp_v_dc*w_vdc*zeta_dc; 
Pref = 1;

Id_limit = 2*sqrt(3/2);

Lf = 0.05;

% ac-link
w_tvc = 0.2*2*pi;
m_v = 0;
Vac_ref = cos(asin(2*Xg*Pref/Ug^2)/2)*Ug;

%system
global system;  
system = "DVCslow"; %"DVCslow" "PLLslow" "overdamp" "DVCslow_TVC"  "PLLslow_TVC" "overdamp_TVC" "TVCfast" "PLLfast_3D" "PLL"

global fault_type; %line_cut voltage_sag frequency
fault_type = "voltage_sag";%"line_cut";"voltage_sag""power_jump""phase_jump"

global model; %line_cut voltage_sag frequency
model = "precise";%"rough";%"precise";%"precise2" "normal" "rough""precise"

global IsTVCslow; %line_cut voltage_sag frequency
IsTVCslow = 0; %0: Constant Iq  1: fault clear Iq

global R_model; 
R_model = "precise";%"num_kp";%"precise";%"normal";%"num1" "normal" "precise"

global R_model_TVC; 
R_model_TVC = "normal2";%"num1" "normal" "precise"


switch fault_type
    case "voltage_sag"
        %fault sag
        Ug_fault = 0.9;  %0.9
        t_c = 0.15 ;    %0.2
    case "line_cut"
    %fault line cut 
        t_c = 0.2;
        X1 = 0.1;
        R1 = 0.01;
        Xgg = Xg-X1;
        Rgg = Rg - R1;
        Xg0 = Xgg/2+X1;
        Rg0 = Rgg/2+R1;
        position = 1; %fault to inf bus
        Rf = 1e-5/(690^2/15e6);
        Im_temp = Rf*(Xgg*position*1j+Rgg*position)/(Rf + Xgg*position*1j+Rgg*position)+Xgg*(1-position)*1j+Rgg*(1-position);
        Imgf = Im_temp *(Xgg*1j+Rgg)/(Xgg*1j+Rgg+Im_temp); %(Rf//Xg*location+Xg*(1-location))//Xg
        Xg_f = imag(Imgf)+X1;
        Lg_f = Xg_f/Ws;
        Rg_f = real(Imgf)+R1;
        Iop_tmp = 1/( ((2-position)*(Xgg*1j+Rgg)*(position)*(Xgg*1j+Rgg))/(2*Xgg*1j+2*Rgg) + Rf); %Xg*(2-location)//Xg*(location)+Rf
        Ug_fault = Iop_tmp*Rf+Iop_tmp/2*position*(1-position)*(Xgg*1j+Rgg);
        Ug_fault_angle = angle(Ug_fault);
        Ug_fault= abs(Ug_fault);
    case "power_jump"
        Pref0 = 0.45;
        t_c = 0.01;
    case "phase_jump"
        phasejump = -18*pi/180;
        t_c = 0.01;
end


%% 
try
    system;
catch
    system = "DVCslow";
end


%% full-order model
%% SEP 
clear xep flag v V Lambda A sig m
x= (0:0.1:1)*2*pi; %delta
y= (0:0.1:1)*5; %Id
n = length(x);
x_set = zeros(5,n^2);
for m = 0:(n^2 - 1)
    index = floor(m);
    index = mod(index,n)+1;
    x_set(1,m+1) = x(index);
    index = floor(m/n);
    index = mod(index,n)+1;
    x_set(4,m+1) = y(index);
    x_set(5,m+1) = 0;
end
torralence = 1e-2; 
mm = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',100000,'OptimalityTolerance',1e-10,'Display','off');
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f_full_pre,xep,options);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set(mm).xep = xep; %#ok<*SAGROW> 
            ep_set(mm).A = A;
            ep_set(mm).Lambda = Lambda;
            ep_set(mm).V = V;   
            ep_set(mm).v = v;     % stable eigenvectors of unstable ep 
            ep_set(mm).flag = flag;
           
            mm = mm+1;
            if flag == 0
                jacob=A;
                prefault_SEP=xep;
            end
        end
    end
end

for mm = 1:length(ep_set)
    disp_v('Index',mm);
    disp_v('Equilibrium',ep_set(mm).xep);
    disp_v('Eigenvalue', ep_set(mm).Lambda);
    disp_v('Eigenvector',ep_set(mm).V);
end



%% reduced order model
if IsTVCslow == 1
   Iq = Iq_pre(end);   % Iq is set as steady state
else
   Iq = Iq00;
end
clear xep flag v V Lambda A sig mm
switch system
    case {"DVCslow","DVCslow_TVC"}
        x=0:0.02:5;
        n = length(x);
        x_set = zeros(2,n);
        x_set(1,:) = x;
    case {"PLLslow", "PLLslow_TVC"}
        x=(-1:0.05:1)*pi/2;
        n = length(x);
        x_set = zeros(2,n);
        x_set(1,:) = x;
    case "overdamp"
        x= (0:0.02:1)*2*pi; %delta
        y= (-0.5:0.02:1)*8; %Id
        n = length(x);
        x_set = zeros(2,n^2);
        for m = 0:(n^2 - 1)
            index = floor(m);
            index = mod(index,n)+1;
            x_set(1,m+1) = x(index);
            index = floor(m/n);
            index = mod(index,n)+1;
            x_set(2,m+1) = y(index);
        end
    case "PLLfast_3D"
            x=0:0.02:5;
            n = length(x);
            x_set = zeros(3,n);
            x_set(1,:) = x;
            m = 1;
            ep_set = [];
            options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',10000,'OptimalityTolerance',1e-10);
            for n = 1:length(x_set(1,:))
                xep = x_set(:,n);
                [xep,ferr,~,~,A] = fsolve(@f,xep,options);
            
                %xep = mod(xep,2*pi);
                torralence = 1e-2; 
                if maxabs(ferr) < torralence
                    if isnewxep(ep_set,xep,torralence)
                       
                        [V,Lambda]=eig(A);
                        Lambda = diag(Lambda);
                        sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
                        sig = (sig + 1)/2;                  % [0,1]
                        flag = sum(sig);                    % number of non-negative eigenvalues
            
                        v = V(:,~sig);                      % the stable sub-space
                        
                        ep_set(m).xep = real(xep); %#ok<*SAGROW> 
                        ep_set(m).A = A;
                        ep_set(m).Lambda = Lambda;
                        ep_set(m).V = V;   
                        ep_set(m).v = v;     % stable eigenvectors of unstable ep 
                        ep_set(m).flag = flag;
                       
                        m = m+1;
                        if flag == 0
                        jacob=A;
                        xeps=xep;
                        end
                    end
                end
            end
            
            for m = 1:length(ep_set)
                disp_v('Index',m);
                disp_v('Equilibrium',ep_set(m).xep);
                disp_v('Eigenvalue', ep_set(m).Lambda);
                disp_v('Eigenvector',ep_set(m).V);
            end
            for n = 1:length(ep_set)
                ep_set_ext(n)=ep_set(n); %#ok<*AGROW> 
            end
            figure;
            color_code = {'blue','magenta','green','red','cyan','blue','magenta','green','red','cyan'};
            color_code2 = {[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0.8500 0.3250 0.0980],'red','cyan','blue','magenta','green','red','cyan'};
            for m = 1:length(ep_set_ext)
                xep = ep_set_ext(m).xep;
                flag = ep_set_ext(m).flag;
                scatter3(xep(1),xep(2),xep(3),color_code{flag+1});
                hold on;
            end    
            n = 1;
            for m = 1 : length(ep_set_ext)        
                flag = ep_set_ext(m).flag;
                if flag 
                    xep = ep_set_ext(m).xep;
                    v = ep_set_ext(m).v;
                    perturb = 5e-1;
            
                    if flag == 1
                        for alpha = (0:0.005:1)*2*pi
                            vp = v(:,1)*sin(alpha) + v(:,2)*cos(alpha);
                            [~ , x_all] = ode45(@f_backward,[0,1],xep+vp*perturb);
                            plot3(x_all(:,1),x_all(:,2),x_all(:,3),'color',color_code2{n});
                        end
                        n = n + 1;
                    elseif flag == 2
                        for beta = [-1,1]
                            vp = v*beta;
                            [~ , x_all] = ode45(@f_backward,[0,0.5],xep+vp*perturb);
                            plot3(x_all(:,1),x_all(:,2),x_all(:,3),'color',color_code2{n});
                        end
                        n = n + 1;
                    end
                    
                end
            end
            axis([0 2.5 -5 5 -4 4]);
end
%%
%try
torralence = 1e-2; 
mm = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',100000,'OptimalityTolerance',1e-10,'Display','off');
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set(mm).xep = xep; %#ok<*SAGROW> 
            ep_set(mm).A = A;
            ep_set(mm).Lambda = Lambda;
            ep_set(mm).V = V;   
            ep_set(mm).v = v;     % stable eigenvectors of unstable ep 
            ep_set(mm).flag = flag;
           
            mm = mm+1;
        end
    end
end

%catch
%end
for mm = 1:length(ep_set)
    disp_v('Index',mm);
    disp_v('Equilibrium',ep_set(mm).xep);
    disp_v('Eigenvalue', ep_set(mm).Lambda);
    disp_v('Eigenvector',ep_set(mm).V);
end

clear ep_set_ext;
switch system
    case {"DVCslow","DVCslow_TVC"}
            ep_set_ext = ep_set; 
    case {"PLLslow","PLLslow_TVC"}
        for n = 1:length(ep_set)
            mm = (n-1)*3;
            ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
            ep_set_ext(mm+2)=ep_set(n);
            ep_set_ext(mm+3)=ep_set(n); 
            ep_set_ext(mm+2).xep(1) = ep_set(n).xep(1) - 2*pi;
            ep_set_ext(mm+3).xep(1) = ep_set(n).xep(1) + 2*pi;
        end
    case "overdamp"
        for n = 1:length(ep_set)
            mm = (n-1)*3;
            ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
            ep_set_ext(mm+2)=ep_set(n);
            ep_set_ext(mm+3)=ep_set(n); 
            ep_set_ext(mm+2).xep(1) = ep_set(n).xep(1) - 2*pi;
            ep_set_ext(mm+3).xep(1) = ep_set(n).xep(1) + 2*pi;
        end
end

f1 = figure(1);
hold on;
grid on;
color_code = {'blue','magenta','red','black'};
switch system
    case {"DVCslow","DVCslow_TVC"}
        ymin= -1.5;%-1.5;
        ymax= 3.5;
        axis([0,2.5,ymin,ymax]);%axis([0,2.5,ymin,ymax]);
        set(gca, 'FontSize', 16);
        for mm = 1 : length(ep_set_ext)
            xep = ep_set_ext(mm).xep;
            flag= ep_set_ext(mm).flag;
            plot(xep(1),xep(2),'o','color',color_code{flag+1},'MarkerSize',8,'linewidth',2);
            %scatter(xep(1),xep(2),color_code{flag+1},'linewidth',2);
            if flag == 1
                v = ep_set_ext(mm).v;
                perturb = 1e-3;
                [~ , x_p] = ode78(@f_backward,[0,30],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,30],xep-v*perturb,odeset('RelTol',1e-5));
                x_all = [flip(x_n,1);x_p];
                plot(x_all(:,1),x_all(:,2),'k-','linewidth',2);
            end
        end
    case {"PLLslow","PLLslow_TVC"}
        ymin=-0.025;
        ymax=0.205;
        axis([0,pi/2,ymin,ymax]);
        xticks(-pi/2:pi/8:pi/2);
        %xticklabels({'$-\frac{\pi}{4}$', '','$0$', '','$\frac{\pi}{4}$','','$\frac{\pi}{2}$'});
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gca, 'FontSize', 16);
        for mm = 1 : length(ep_set_ext)
            xep = ep_set_ext(mm).xep;
            flag= ep_set_ext(mm).flag;
            plot(xep(1),xep(2),'o','color',color_code{flag+1},'MarkerSize',8,'linewidth',2);
            %scatter(xep(1),xep(2),color_code{flag+1},'linewidth',2);
            if flag == 1
                v = ep_set_ext(mm).v;
                perturb = 1e-3;
                [~ , x_p] = ode78(@f_backward,[0,25],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,25],xep-v*perturb,odeset('RelTol',1e-5));
                x_all = [flip(x_n,1);x_p];
                plot(x_all(:,1),x_all(:,2),'k-','linewidth',2);
            end
        end
    case "overdamp"
        ymin=-2;
        ymax=8;
        axis([-2*pi,2*pi,ymin,ymax]);
        xticks(-2*pi:pi/2:2*pi);
        xticklabels({'$-2\pi$','','$-\pi$', '','$0$', '','$\pi$','','$2\pi$'});
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gca, 'FontSize', 14);
        for mm = 1 : length(ep_set_ext)
            xep = ep_set_ext(mm).xep;
            flag= ep_set_ext(mm).flag;
            scatter(xep(1),xep(2),color_code{flag+1});
            if flag == 1
                v = ep_set_ext(mm).v;
                perturb = 1e-3;
                [~ , x_p] = ode78(@f_backward,[0,20],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,20],xep-v*perturb,odeset('RelTol',1e-5));
                x_all = [flip(x_n,1);x_p];
                plot(x_all(:,1),x_all(:,2),'k-','linewidth',2);
            end
        end
end

%%  fault clear time DOA  ( if considering slow TVC run this part)
if IsTVCslow == 1
Iq = Iq_clear;  % Iq is set as clearing point value
clear xep flag v V Lambda A sig mm
ep_set = [];
switch system
    case "DVCslow"
        x=0:0.02:5;
        n = length(x);
        x_set = zeros(2,n);
        x_set(1,:) = x;
    case"PLLslow"
        x=(-1:0.05:1)*pi/2;
        n = length(x);
        x_set = zeros(2,n);
        x_set(1,:) = x;
    case "overdamp"
        x= (0:0.02:1)*2*pi; %delta
        y= (-0.5:0.02:1)*8; %Id
        n = length(x);
        x_set = zeros(2,n^2);
        for m = 0:(n^2 - 1)
            index = floor(m);
            index = mod(index,n)+1;
            x_set(1,m+1) = x(index);
            index = floor(m/n);
            index = mod(index,n)+1;
            x_set(2,m+1) = y(index);
        end
end
mm = 1;
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set(mm).xep = xep; %#ok<*SAGROW> 
            ep_set(mm).A = A;
            ep_set(mm).Lambda = Lambda;
            ep_set(mm).V = V;   
            ep_set(mm).v = v;     % stable eigenvectors of unstable ep 
            ep_set(mm).flag = flag;
           
            mm = mm+1;
        end
    end
end

clear ep_set_ext;
switch system
    case "DVCslow"
            ep_set_ext = ep_set; 
    case "PLLslow"
        for n = 1:length(ep_set)
            mm = (n-1)*3;
            ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
            ep_set_ext(mm+2)=ep_set(n);
            ep_set_ext(mm+3)=ep_set(n); 
            ep_set_ext(mm+2).xep(1) = ep_set(n).xep(1) - 2*pi;
            ep_set_ext(mm+3).xep(1) = ep_set(n).xep(1) + 2*pi;
        end
    case "overdamp"
        for n = 1:length(ep_set)
            mm = (n-1)*3;
            ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
            ep_set_ext(mm+2)=ep_set(n);
            ep_set_ext(mm+3)=ep_set(n); 
            ep_set_ext(mm+2).xep(1) = ep_set(n).xep(1) - 2*pi;
            ep_set_ext(mm+3).xep(1) = ep_set(n).xep(1) + 2*pi;
        end
end

figure(f1);
hold on;
switch system
    case "DVCslow"
        ymin= -1.5;
        ymax= 3;
        axis([0,2.5,ymin,ymax]);
        for mm = 1 : length(ep_set_ext)
            xep = ep_set_ext(mm).xep;
            flag= ep_set_ext(mm).flag;
            scatter(xep(1),xep(2),color_code{flag+1});
            if flag == 1
                v = ep_set_ext(mm).v;
                perturb = 1e-3;
                [~ , x_p] = ode78(@f_backward,[0,20],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,20],xep-v*perturb,odeset('RelTol',1e-5));
                x_all = [flip(x_n,1);x_p];
                plot(x_all(:,1),x_all(:,2),'m-','linewidth',2);
            end
        end
    case "PLLslow"
        ymin=-1;
        ymax=1;
        axis([-pi/2,pi/2,ymin,ymax]);
        xticks(-pi/2:pi/4:pi/2);
        xticklabels({'$-\frac{\pi}{2}$', '','$0$', '','$\frac{\pi}{2}$'});
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gca, 'FontSize', 14);
        for mm = 1 : length(ep_set_ext)
            xep = ep_set_ext(mm).xep;
            flag= ep_set_ext(mm).flag;
            scatter(xep(1),xep(2),color_code{flag+1});
            if flag == 1
                v = ep_set_ext(mm).v;
                perturb = 1e-3;
                [~ , x_p] = ode78(@f_backward,[0,20],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,20],xep-v*perturb,odeset('RelTol',1e-5));
                x_all = [flip(x_n,1);x_p];
                plot(x_all(:,1),x_all(:,2),'m-','linewidth',2);
            end
        end
    case "overdamp"
        ymin=-2;
        ymax=8;
        axis([-2*pi,2*pi,ymin,ymax]);
        xticks(-2*pi:pi/2:2*pi);
        xticklabels({'$-2\pi$','','$-\pi$', '','$0$', '','$\pi$','','$2\pi$'});
        set(gca, 'TickLabelInterpreter', 'latex');
        set(gca, 'FontSize', 14);
        for mm = 1 : length(ep_set_ext)
            xep = ep_set_ext(mm).xep;
            flag= ep_set_ext(mm).flag;
            scatter(xep(1),xep(2),color_code{flag+1});
            if flag == 1
                v = ep_set_ext(mm).v;
                perturb = 1e-3;
                [~ , x_p] = ode78(@f_backward,[0,20],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,20],xep-v*perturb,odeset('RelTol',1e-5));
                x_all = [flip(x_n,1);x_p];
                plot(x_all(:,1),x_all(:,2),'m-','linewidth',2);
            end
        end
end

end

%% time domain   %% order-reduced model results according to fault type
switch fault_type
    %% voltage sag
    case "voltage_sag"
    t_start = 0.1;
    delta_pre = [prefault_SEP(1); prefault_SEP(1)];
    Int_pre = [prefault_SEP(2); prefault_SEP(2)];
    ydc_pre = [prefault_SEP(3); prefault_SEP(3)];
    Intdc_pre = [prefault_SEP(4); prefault_SEP(4)];
    Iq_pre = [prefault_SEP(5); prefault_SEP(5)];
    Vdc_pre = sqrt(ydc_pre+Vdc_ref^2);
    Id_pre = Intdc_pre + kp_v_dc*ydc_pre;
    Vq_pre = (Xg*Id_pre+Rg*Iq-Ug*sin(delta_pre)+Id_pre.*Lg.*Int_pre)./(1-Id_pre*Lg*kp_pll);
    omega_pre=kp_pll.*Vq_pre+Int_pre;
    t_prefault = [0;t_start];


    [t_fault , x_all] = ode78(@f_fault,[t_start,t_start+t_c],prefault_SEP,odeset('RelTol',1e-10));
    delta_fault = x_all(:,1);
    Int_fault = x_all(:,2);
    ydc_fault = x_all(:,3);
    Intdc_fault = x_all(:,4);
    Iq_fault = x_all(:,5);
    Vdc_fault = sqrt(ydc_fault+Vdc_ref^2);
    Id_fault = Intdc_fault + kp_v_dc*ydc_fault;
    Vq_fault = (Xg*Id_fault+Rg*Iq_fault-Ug_fault*sin(delta_fault)+Id_fault.*Lg.*Int_fault)./(1-Id_fault*Lg*kp_pll);
    omega_fault=kp_pll.*Vq_fault+Int_fault;

  
    Iq_clear = Iq_fault(end);


    [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,:),odeset('RelTol',1e-10));
    delta_post = x_all2(:,1);
    Int_post = x_all2(:,2);
    ydc_post = x_all2(:,3);
    Intdc_post = x_all2(:,4);
    Iq_post = x_all2(:,5);
    Vdc_post = sqrt(ydc_post+Vdc_ref^2);
    Id_post = Intdc_post + kp_v_dc*ydc_post;
    Vq_post = (Xg*Id_post+Rg*Iq_post-Ug*sin(delta_post)+Id_post.*Lg.*Int_post)./(1-Id_post*Lg*kp_pll);
    omega_post=kp_pll.*Vq_post+Int_post;

    % fault-clear t0
    [t_fault_clear , x_all3] = ode78(@f_post,[t_fault(end),t_fault(end)+t_0],x_all(end,:),odeset('RelTol',1e-10));
    delta_fault_clear = x_all3(:,1);
    Int_fault_clear = x_all3(:,2);
    ydc_fault_clear = x_all3(:,3);
    Intdc_fault_clear = x_all3(:,4);
    Iq_fault_clear = x_all3(:,5);
    Vdc_fault_clear = sqrt(ydc_fault_clear+Vdc_ref^2);
    Id_fault_clear = Intdc_fault_clear + kp_v_dc*ydc_fault_clear;
    Vq_fault_clear = (Xg*Id_fault_clear+Rg*Iq_fault_clear-Ug*sin(delta_fault_clear)+Id_fault_clear.*Lg.*Int_fault_clear)./(1-Id_fault_clear*Lg*kp_pll);
    omega_fault_clear=kp_pll.*Vq_fault_clear+Int_fault_clear;
    
    % reduce order
    if IsTVCslow == 1
        Iq = Iq_clear;  %for slow TVC
    else
        Iq = Iq00;
    end
    switch system
        case "DVCslow"
        % during fault RM2
        temp_ini = [delta_pre;0];
        delta_fault_RM2 = asin((Id_fault*Xg+Iq*Rg)/Ug_fault);

        %postfault RM
        Intdc0 = Intdc_fault_clear(end);
        ydc0 = ydc_fault_clear(end);

        %boundary layer correction
        Id_s0 = Intdc0 + kp_v_dc*ydc0;
        delta_s0 = asin(Id_s0*Xg/Ug);
        delta_0 = delta_fault_clear(end);
        ydc_star = 2/C_dc*Id_s0*sin(delta_s0)*(delta_0-delta_s0)/cos(delta_s0)/kp_pll*(1-kp_pll*Lg*Id_s0);

        [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_fault(end)+t_0:0.001:t_end],[ydc0+ydc_star;Intdc0],odeset('RelTol',1e-10));
        Intdc_post_RM = x_all3(:,2);
        ydc_post_RM = x_all3(:,1);
        Id_post_RM = Intdc_post_RM + kp_v_dc*ydc_post_RM;
        Vdc_post_RM = sqrt(ydc_post_RM+Vdc_ref^2);
        delta_post_RM = asin((Xg*Id_post_RM+Rg*Iq)/Ug);
        
        % considering t0
        Intdc_post_RMc = Intdc_post_RM; %[Intdc_fault_clear; Intdc_post_RM ];
        ydc_post_RMc = ydc_post_RM;%[ydc_fault_clear; ydc_post_RM];
        Id_post_RMc = Intdc_post_RMc + kp_v_dc*ydc_post_RMc;
        delta_post_RMc = asin((Xg*Id_post_RMc+Rg*Iq)/Ug);
        Vdc_post_RMc = sqrt(ydc_post_RMc+Vdc_ref^2);
        t_postfault_RMc = t_postfault_RM;%[t_fault_clear; t_postfault_RM];


        case "PLLslow"
        Id_fault_RM2 = Pref./Ug_fault./cos(delta_fault);

        %postfault RM
        Int0 = Int_fault_clear(end);
        delta0 = delta_fault_clear(end);

        Id_slow0 = (Pref+Iq.*Ug.*sin(delta0)-Iq.^2*Rg)./Ug./cos(delta0) - (Pref+Iq.*Ug.*sin(delta0)-Iq.^2*Rg).^2./Ug^3./cos(delta0).^3*Rg;
        ww = sqrt(Ug*cos(delta0)-Ug^2*cos(delta0)^2)/2;
        C1 = Id_fault_clear(end)-Id_slow0;
        C2 = (-C1*Ug*cos(delta0)/2+C_dc/2*zeta_dc*w_vdc*ydc_fault_clear(end))/ww;
        Id_star = (C1*Ug*cos(delta0)/2+C2*ww)/(Ug^2*cos(delta0)^2/4+ww^2);
        Id_star=Id_star*1.3;
        delta_star = 1/w_vdc*Id_star*(Xg*kp_pll)/(1-kp_pll*Lg*Id_slow0);
        Int_star = 1/w_vdc*Id_star*(Xg*ki_pll)/(1-kp_pll*Lg*Id_slow0);


        [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_fault(end)+t_0,t_end],[delta0+delta_star;Int0+Int_star],odeset('RelTol',1e-10));
        Int_post_RM = x_all3(:,2);
        delta_post_RM = x_all3(:,1);
        Id_post_RM = (Pref+Iq.*Ug.*sin(delta_post_RM)-Iq.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq.*Ug.*sin(delta_post_RM)-Iq.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;
    

        % considering t0
        Int_post_RMc = Int_post_RM;
        delta_post_RMc = delta_post_RM;
        Id_post_RMc = (Pref+Iq.*Ug.*sin(delta_post_RMc)-Iq.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq.*Ug.*sin(delta_post_RMc)-Iq.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
        t_postfault_RMc = t_postfault_RM;

        case {"TVCfast", "PLLslow_TVC", "DVCslow_TVC"}
        clear x_all x_all2
    
        [t_fault_Rtvc , x_all] = ode78(@f_TVCfast_fault_time,[t_start,t_start+t_c],prefault_SEP(1:4),odeset('RelTol',1e-10));
        delta_fault_Rtvc = x_all(:,1);
        Int_fault_Rtvc = x_all(:,2);
        ydc_fault_Rtvc = x_all(:,3);
        Intdc_fault_Rtvc = x_all(:,4);   
        Vdc_fault_Rtvc = sqrt(ydc_fault_Rtvc+Vdc_ref^2);
        Id_fault_Rtvc = Intdc_fault_Rtvc + kp_v_dc*ydc_fault_Rtvc;
        Iq_fault_Rtvc = (Ug_fault*cos(delta_fault_Rtvc)+Id_fault_Rtvc*Rg-Vac_ref)/(Xg+1/m_v);
        Vq_fault_Rtvc = (Xg*Id_fault_Rtvc+Rg*Iq_fault_Rtvc-Ug_fault*sin(delta_fault_Rtvc)+Id_fault_Rtvc.*Lg.*Int_fault_Rtvc)./(1-Id_fault_Rtvc*Lg*kp_pll);
        omega_fault_Rtvc=kp_pll.*Vq_fault_Rtvc+Int_fault_Rtvc;
        
        [t_postfault_Rtvc , x_all2] = ode78(@f_TVCfast_post_time,[t_fault_Rtvc(end),t_end],x_all(end,:),odeset('RelTol',1e-10));
        delta_post_Rtvc = x_all2(:,1);
        Int_post_Rtvc = x_all2(:,2);
        ydc_post_Rtvc = x_all2(:,3);
        Intdc_post_Rtvc = x_all2(:,4);
        Vdc_post_Rtvc = sqrt(ydc_post_Rtvc+Vdc_ref^2);
        Id_post_Rtvc = Intdc_post_Rtvc + kp_v_dc*ydc_post_Rtvc;
        Iq_post_Rtvc = (Ug*cos(delta_post_Rtvc)+Id_post_Rtvc*Rg-Vac_ref)/(Xg+1/m_v);
        Vq_post_Rtvc = (Xg*Id_post_Rtvc+Rg*Iq_post_Rtvc-Ug*sin(delta_post_Rtvc)+Id_post_Rtvc.*Lg.*Int_post_Rtvc)./(1-Id_post_Rtvc*Lg*kp_pll);
        omega_post_Rtvc=kp_pll.*Vq_post_Rtvc+Int_post_Rtvc;
    
        t_timedomain_Rtvc = [t_prefault;t_fault_Rtvc;t_postfault_Rtvc];
        delta_timedomain_Rtvc = [delta_pre; delta_fault_Rtvc; delta_post_Rtvc];
        Int_timedomain_Rtvc = [Int_pre; Int_fault_Rtvc; Int_post_Rtvc];%./2/pi+50;
        Vdc_timedomain_Rtvc = [Vdc_pre; Vdc_fault_Rtvc; Vdc_post_Rtvc];
        Id_timedomain_Rtvc = [Id_pre; Id_fault_Rtvc; Id_post_Rtvc];
        Iq_timedomain_Rtvc = [Iq_pre; Iq_fault_Rtvc; Iq_post_Rtvc];
        omega_timedomain_Rtvc = [omega_pre; omega_fault_Rtvc; omega_post_Rtvc];

        % fault-clear t0
        [t_fault_clear_Rtvc , x_all3] = ode78(@f_TVCfast_post_time,[t_fault_Rtvc(end),t_fault_Rtvc(end)+t_0],x_all(end,:),odeset('RelTol',1e-10));
        delta_fault_clear_Rtvc = x_all3(:,1);
        Int_fault_clear_Rtvc = x_all3(:,2);
        ydc_fault_clear_Rtvc = x_all3(:,3);
        Intdc_fault_clear_Rtvc = x_all3(:,4);
        Vdc_fault_clear_Rtvc = sqrt(ydc_fault_clear_Rtvc+Vdc_ref^2);
        Id_fault_clear_Rtvc = Intdc_fault_clear_Rtvc + kp_v_dc*ydc_fault_clear_Rtvc;
        Iq_fault_clear_Rtvc = (Ug*cos(delta_fault_clear_Rtvc)+Id_fault_clear_Rtvc*Rg-Vac_ref)/(Xg+1/m_v);
        Vq_fault_clear_Rtvc = (Xg*Id_fault_clear_Rtvc+Rg*Iq_fault_clear_Rtvc-Ug*sin(delta_fault_clear_Rtvc)+Id_fault_clear_Rtvc.*Lg.*Int_fault_clear_Rtvc)./(1-Id_fault_clear_Rtvc*Lg*kp_pll);
        omega_fault_clear_Rtvc=kp_pll.*Vq_fault_clear_Rtvc+Int_fault_clear_Rtvc;
        


        if system == "DVCslow_TVC"
            % during fault RM2
            Intdc0 = prefault_SEP(4);
            ydc0 = prefault_SEP(3);
            [t_postfault_RM2, x_all3] = ode78(@f_reduce_fault,[t_start,t_start+t_c],[ydc0;Intdc0],odeset('RelTol',1e-10));
            Intdc_fault_RM2 = x_all3(:,2);
            ydc_fault_RM2 = x_all3(:,1);
            Id_fault_RM2 = Intdc_fault_RM2 + kp_v_dc*ydc_fault_RM2;
            Vdc_fault_RM2 = sqrt(ydc_fault_RM2+Vdc_ref^2);
            delta_fault_RM2 = asin((Xg*Id_fault_RM2)/Ug_fault);
            Iq_fault_RM2 = (Ug_fault*cos(delta_fault_RM2) + Id_fault_RM2*Rg - Vac_ref)*m_v/(Xg*m_v+1);
            delta_fault_RM2 = asin((Id_fault_RM2*Xg+Iq_fault_RM2*Rg)/Ug_fault);
            Iq_fault_RM2 = (Ug_fault*cos(delta_fault_RM2) + Id_fault_RM2*Rg - Vac_ref)*m_v/(Xg*m_v+1);

            delta_fault_RM22 = asin((Id_fault*Xg)/Ug_fault);
            Iq_fault_RM2 = (Ug_fault*cos(delta_fault_RM22) + Id_fault*Rg - Vac_ref)*m_v/(Xg*m_v+1);
            delta_fault_RM22 = asin((Id_fault*Xg+Iq_fault_RM2*Rg)/Ug_fault);
            Iq_fault_RM2 = (Ug_fault*cos(delta_fault_RM22) + Id_fault*Rg - Vac_ref)*m_v/(Xg*m_v+1);



            
            % postfault RM
            Intdc0 = Intdc_fault_clear(end);
            ydc0 = ydc_fault_clear(end);

            % boundary layer correction
            Id_s0 = Intdc0 + kp_v_dc*ydc0;
            delta_s0 = asin(Id_s0*Xg/Ug);
            delta_0 = delta_fault_clear(end);
            ydc_star = 2/C_dc*Id_s0*sin(delta_s0)*(delta_0-delta_s0)/cos(delta_s0)/kp_pll*(1-kp_pll*Lg*Id_s0);

            [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_fault_Rtvc(end)+t_0,t_end],[ydc0+ydc_star;Intdc0],odeset('RelTol',1e-10));
            Intdc_post_RM = x_all3(:,2);
            ydc_post_RM = x_all3(:,1);
            Id_post_RM = Intdc_post_RM + kp_v_dc*ydc_post_RM;
            Vdc_post_RM = sqrt(ydc_post_RM+Vdc_ref^2);
            delta_post_RM = asin((Xg*Id_post_RM)/Ug);
            Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg - Vac_ref)*m_v/(Xg*m_v+1);
            delta_post_RM = asin((Id_post_RM*Xg+Iq_post_RM*Rg)/Ug);
            Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg - Vac_ref)*m_v/(Xg*m_v+1);


            % considering t0
            Intdc_post_RMc = [Intdc_post_RM ];
            ydc_post_RMc = [ydc_post_RM];
            Id_post_RMc = Intdc_post_RMc + kp_v_dc*ydc_post_RMc;
            Vdc_post_RMc = sqrt(ydc_post_RMc+Vdc_ref^2);
            delta_post_RMc = asin((Xg*Id_post_RMc)/Ug);
            Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg - Vac_ref)*m_v/(Xg*m_v+1);
            delta_post_RMc = asin((Id_post_RMc*Xg+Iq_post_RMc*Rg)/Ug);
            Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg - Vac_ref)*m_v/(Xg*m_v+1);
            t_postfault_RMc = [ t_postfault_RM];

        end
        if system == "PLLslow_TVC"
            Iq_fault_RM2 = (Ug_fault*cos(delta_fault) - Vac_ref)*m_v/(1+Xg*m_v);
            Id_fault_RM2 = (Pref+Iq_fault_RM2.*Ug_fault.*sin(delta_fault)-Iq_fault_RM2.^2*Rg)./Ug_fault./cos(delta_fault) - (Pref+Iq_fault_RM2.*Ug_fault.*sin(delta_fault)-Iq_fault_RM2.^2*Rg).^2./Ug_fault^3./cos(delta_fault).^3*Rg;
            Iq_fault_RM2 = (Ug_fault*cos(delta_fault) + Id_fault_RM2*Rg- Vac_ref)*m_v/(1+Xg*m_v);
            Id_fault_RM2 = (Pref+Iq_fault_RM2.*Ug_fault.*sin(delta_fault)-Iq_fault_RM2.^2*Rg)./Ug_fault./cos(delta_fault) - (Pref+Iq_fault_RM2.*Ug_fault.*sin(delta_fault)-Iq_fault_RM2.^2*Rg).^2./Ug_fault^3./cos(delta_fault).^3*Rg;

            %boundary layer correction
            delta0 = delta_fault_clear(end);
            Iqq = (Ug*cos(delta0) - Vac_ref )*m_v/(1+Xg*m_v);
            Id_slow0 = (Pref+Iqq.*Ug.*sin(delta0))./Ug./cos(delta0);
            ww = sqrt(Ug*cos(delta0)-Ug^2*cos(delta0)^2)/2;
            C1 = Id_fault_clear(end)-Id_slow0;
            C2 = (-C1*Ug*cos(delta0)/2+C_dc/2*zeta_dc*w_vdc*ydc_fault_clear(end))/ww;
            Id_star = (C1*Ug*cos(delta0)/2+C2*ww)/(Ug^2*cos(delta0)^2/4+ww^2);
            Id_star=Id_star*1.3;
            delta_star = 1/w_vdc*Id_star*(Xg*kp_pll)/(1-kp_pll*Lg*Id_slow0);
            Int_star = 1/w_vdc*Id_star*(Xg*ki_pll)/(1-kp_pll*Lg*Id_slow0);

            %postfault RM
            Int0 = Int_fault_clear_Rtvc(end);
            delta0 = delta_fault_clear_Rtvc(end);
            [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_fault(end)+t_0,t_end],[delta0+delta_star;Int0+Int_star],odeset('RelTol',1e-10));
            Int_post_RM = x_all3(:,2);
            delta_post_RM = x_all3(:,1);
            Iq_post_RM = (Ug*cos(delta_post_RM) - Vac_ref)*m_v/(1+Xg*m_v);
            Id_post_RM = (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;
            Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg- Vac_ref)*m_v/(1+Xg*m_v);
            Id_post_RM = (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;

            % considering t0
            Int_post_RMc = [Int_fault_clear_Rtvc; Int_post_RM ];
            delta_post_RMc = [delta_fault_clear_Rtvc; delta_post_RM];
            Iq_post_RMc = (Ug*cos(delta_post_RMc) - Vac_ref)*m_v/(1+Xg*m_v);
            Id_post_RMc = (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
            Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg- Vac_ref)*m_v/(1+Xg*m_v);
            Id_post_RMc = (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
            t_postfault_RMc = [t_fault_clear_Rtvc; t_postfault_RM];

        end


    end
    %% line cut
    case "line_cut"
    

    %% power jump
    case "power_jump"
        t_start = 0.1;
        delta_pre = [prefault_SEP(1); prefault_SEP(1)];
        Int_pre = [prefault_SEP(2); prefault_SEP(2)];
        ydc_pre = [prefault_SEP(3); prefault_SEP(3)];
        Intdc_pre = [prefault_SEP(4); prefault_SEP(4)];
        Iq_pre = [prefault_SEP(5); prefault_SEP(5)];
        Vdc_pre = sqrt(ydc_pre+Vdc_ref^2);
        Id_pre = Intdc_pre + kp_v_dc*ydc_pre;
        Vq_pre = (Xg*Id_pre+Rg*Iq-Ug*sin(delta_pre)+Id_pre.*Lg.*Int_pre)./(1-Id_pre*Lg*kp_pll);
        omega_pre=kp_pll.*Vq_pre+Int_pre;
        t_prefault = [0;t_start];

        Iq_clear = Iq_pre(end);

        t_fault = [];
        delta_fault = [];
        Int_fault = [];
        ydc_fault = [];
        Intdc_fault = [];
        Vdc_fault = [];
        Id_fault = [];
        Vq_fault = [];
        omega_fault=[];
        Iq_fault = [];
      
    
        [t_postfault , x_all2] = ode78(@f_post,[t_start,t_end],prefault_SEP,odeset('RelTol',1e-6));
        delta_post = x_all2(:,1);
        Int_post = x_all2(:,2);
        ydc_post = x_all2(:,3);
        Intdc_post = x_all2(:,4);
        Iq_post = x_all2(:,5);
        Vdc_post = sqrt(ydc_post+Vdc_ref^2);
        Id_post = Intdc_post + kp_v_dc*ydc_post;
        Vq_post = (Xg*Id_post+Rg*Iq-Ug*sin(delta_post)+Id_post.*Lg.*Int_post)./(1-Id_post*Lg*kp_pll);
        omega_post=kp_pll.*Vq_post+Int_post;


        % fault-clear t0
        [t_fault_clear , x_all3] = ode78(@f_post,[t_start,t_start+t_0],prefault_SEP,odeset('RelTol',1e-10));
        delta_fault_clear = x_all3(:,1);
        Int_fault_clear = x_all3(:,2);
        ydc_fault_clear = x_all3(:,3);
        Intdc_fault_clear = x_all3(:,4);
        Iq_fault_clear = x_all3(:,5);
        Vdc_fault_clear = sqrt(ydc_fault_clear+Vdc_ref^2);
        Id_fault_clear = Intdc_fault_clear + kp_v_dc*ydc_fault_clear;
        Vq_fault_clear = (Xg*Id_fault_clear+Rg*Iq_fault_clear-Ug*sin(delta_fault_clear)+Id_fault_clear.*Lg.*Int_fault_clear)./(1-Id_fault_clear*Lg*kp_pll);
        omega_fault_clear=kp_pll.*Vq_fault_clear+Int_fault_clear;
                
        
        % reduce order
        if IsTVCslow == 1
            Iq = Iq_clear;  %for slow TVC
        else
            Iq = Iq00;
        end
        switch system
            case "DVCslow"

                delta_fault_RM2 =[];
                %postfault RM
                Intdc0 = Intdc_fault_clear(end);
                ydc0 = ydc_fault_clear(end);
                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[ydc0;Intdc0],odeset('RelTol',1e-10));
                Intdc_post_RM = x_all3(:,2);
                ydc_post_RM = x_all3(:,1);
                Id_post_RM = Intdc_post_RM + kp_v_dc*ydc_post_RM;
                Vdc_post_RM = sqrt(ydc_post_RM+Vdc_ref^2);
                delta_post_RM = asin((Xg*Id_post_RM+Rg*Iq)/Ug);

                 % considering t0
                Intdc_post_RMc = [Intdc_fault_clear; Intdc_post_RM ];
                ydc_post_RMc = [ydc_fault_clear; ydc_post_RM];
                Id_post_RMc = Intdc_post_RMc + kp_v_dc*ydc_post_RMc;
                delta_post_RMc = asin((Xg*Id_post_RMc+Rg*Iq)/Ug);
                Vdc_post_RMc = sqrt(ydc_post_RMc+Vdc_ref^2);
                t_postfault_RMc = [t_fault_clear; t_postfault_RM];


             case "PLLslow"

                Id_fault_RM2 = [];
                %postfault RM
                Int0 =  Int_fault_clear(end);
                delta0 = delta_fault_clear(end);
                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[delta0;Int0],odeset('RelTol',1e-10));
                Int_post_RM = x_all3(:,2);
                delta_post_RM = x_all3(:,1);
                Id_post_RM = (Pref+Iq.*Ug.*sin(delta_post_RM)-Iq.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq.*Ug.*sin(delta_post_RM)-Iq.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;

                % considering t0
                Int_post_RMc = [Int_fault_clear; Int_post_RM ];
                delta_post_RMc = [delta_fault_clear; delta_post_RM];
                Id_post_RMc = (Pref+Iq.*Ug.*sin(delta_post_RMc)-Iq.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq.*Ug.*sin(delta_post_RMc)-Iq.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
                t_postfault_RMc = [t_fault_clear; t_postfault_RM];



            case {"TVCfast", "PLLslow_TVC", "DVCslow_TVC"}
            clear x_all x_all2
        
           
            delta_fault_Rtvc = [];
            Int_fault_Rtvc = [];
            ydc_fault_Rtvc = [];
            Intdc_fault_Rtvc = [];   
            Vdc_fault_Rtvc = [];
            Id_fault_Rtvc = [];
            Iq_fault_Rtvc = [];
            Vq_fault_Rtvc = [];
            omega_fault_Rtvc=[];
            t_fault_Rtvc=[];
            
            [t_postfault_Rtvc , x_all2] = ode78(@f_TVCfast_post_time,[t_start,t_end],prefault_SEP(1:4),odeset('RelTol',1e-10));
            delta_post_Rtvc = x_all2(:,1);
            Int_post_Rtvc = x_all2(:,2);
            ydc_post_Rtvc = x_all2(:,3);
            Intdc_post_Rtvc = x_all2(:,4);
            Vdc_post_Rtvc = sqrt(ydc_post_Rtvc+Vdc_ref^2);
            Id_post_Rtvc = Intdc_post_Rtvc + kp_v_dc*ydc_post_Rtvc;
            Iq_post_Rtvc = (Ug*cos(delta_post_Rtvc)+Id_post_Rtvc*Rg-Vac_ref)/(Xg+1/m_v);
            Vq_post_Rtvc = (Xg*Id_post_Rtvc+Rg*Iq_post_Rtvc-Ug*sin(delta_post_Rtvc)+Id_post_Rtvc.*Lg.*Int_post_Rtvc)./(1-Id_post_Rtvc*Lg*kp_pll);
            omega_post_Rtvc=kp_pll.*Vq_post_Rtvc+Int_post_Rtvc;
        
            t_timedomain_Rtvc = [t_prefault;t_fault_Rtvc;t_postfault_Rtvc];
            delta_timedomain_Rtvc = [delta_pre; delta_fault_Rtvc; delta_post_Rtvc];
            Int_timedomain_Rtvc = [Int_pre; Int_fault_Rtvc; Int_post_Rtvc];%./2/pi+50;
            Vdc_timedomain_Rtvc = [Vdc_pre; Vdc_fault_Rtvc; Vdc_post_Rtvc];
            Id_timedomain_Rtvc = [Id_pre; Id_fault_Rtvc; Id_post_Rtvc];
            Iq_timedomain_Rtvc = [Iq_pre; Iq_fault_Rtvc; Iq_post_Rtvc];
            omega_timedomain_Rtvc = [omega_pre; omega_fault_Rtvc; omega_post_Rtvc];

            % fault-clear t0
            [t_fault_clear_Rtvc , x_all3] = ode78(@f_TVCfast_post_time,[t_start,t_start+t_0],prefault_SEP(1:4),odeset('RelTol',1e-10));
            delta_fault_clear_Rtvc = x_all3(:,1);
            Int_fault_clear_Rtvc = x_all3(:,2);
            ydc_fault_clear_Rtvc = x_all3(:,3);
            Intdc_fault_clear_Rtvc = x_all3(:,4);
            Vdc_fault_clear_Rtvc = sqrt(ydc_fault_clear_Rtvc+Vdc_ref^2);
            Id_fault_clear_Rtvc = Intdc_fault_clear_Rtvc + kp_v_dc*ydc_fault_clear_Rtvc;
            Iq_fault_clear_Rtvc = (Ug*cos(delta_fault_clear_Rtvc)+Id_fault_clear_Rtvc*Rg-Vac_ref)/(Xg+1/m_v);
            Vq_fault_clear_Rtvc = (Xg*Id_fault_clear_Rtvc+Rg*Iq_fault_clear_Rtvc-Ug*sin(delta_fault_clear_Rtvc)+Id_fault_clear_Rtvc.*Lg.*Int_fault_clear_Rtvc)./(1-Id_fault_clear_Rtvc*Lg*kp_pll);
            omega_fault_clear_Rtvc=kp_pll.*Vq_fault_clear_Rtvc+Int_fault_clear_Rtvc;

            if system == "DVCslow_TVC"      
                %postfault RM
                delta_fault_RM2 =[];
                Iq_fault_RM2 = [];
                Intdc0 = Intdc_pre(end);
                ydc0 = ydc_pre(end);
                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[ydc0;Intdc0],odeset('RelTol',1e-10));
                Intdc_post_RM = x_all3(:,2);
                ydc_post_RM = x_all3(:,1);
                Id_post_RM = Intdc_post_RM + kp_v_dc*ydc_post_RM;
                Vdc_post_RM = sqrt(ydc_post_RM+Vdc_ref^2);
                delta_post_RM = asin((Xg*Id_post_RM)/Ug);
                Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg - Vac_ref)*m_v/(Xg*m_v+1);
                delta_post_RM = asin((Id_post_RM*Xg+Iq_post_RM*Rg)/Ug);
                Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg - Vac_ref)*m_v/(Xg*m_v+1);


                % considering t0
                Intdc_post_RMc = [Intdc_fault_clear_Rtvc; Intdc_post_RM ];
                ydc_post_RMc = [ydc_fault_clear_Rtvc; ydc_post_RM];
                Id_post_RMc = Intdc_post_RMc + kp_v_dc*ydc_post_RMc;
                Vdc_post_RMc = sqrt(ydc_post_RMc+Vdc_ref^2);
                delta_post_RMc = asin((Xg*Id_post_RMc)/Ug);
                Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg - Vac_ref)*m_v/(Xg*m_v+1);
                delta_post_RMc = asin((Id_post_RMc*Xg+Iq_post_RMc*Rg)/Ug);
                Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg - Vac_ref)*m_v/(Xg*m_v+1);
                t_postfault_RMc = [t_fault_clear_Rtvc; t_postfault_RM];

            end
            if system == "PLLslow_TVC"
                Id_fault_RM2 = [];
                Iq_fault_RM2 = [];
                %postfault RM
                delta0 = delta_pre(end);
                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[delta0;Int0],odeset('RelTol',1e-10));
                Int_post_RM = x_all3(:,2);
                delta_post_RM = x_all3(:,1);
                Iq_post_RM = (Ug*cos(delta_post_RM) - Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RM = (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;
                Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg- Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RM = (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;

                % considering t0
                Int_post_RMc = [Int_fault_clear_Rtvc; Int_post_RM ];
                delta_post_RMc = [delta_fault_clear_Rtvc; delta_post_RM];
                Iq_post_RMc = (Ug*cos(delta_post_RMc) - Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RMc = (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
                Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg- Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RMc = (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
                t_postfault_RMc = [t_fault_clear_Rtvc; t_postfault_RM];
            end
        end

    %% phase jump
    case "phase_jump"
        t_start = 0.1;
        delta_pre = [prefault_SEP(1); prefault_SEP(1)];
        Int_pre = [prefault_SEP(2); prefault_SEP(2)];
        ydc_pre = [prefault_SEP(3); prefault_SEP(3)];
        Intdc_pre = [prefault_SEP(4); prefault_SEP(4)];
        Iq_pre = [prefault_SEP(5); prefault_SEP(5)];
        Vdc_pre = sqrt(ydc_pre+Vdc_ref^2);
        Id_pre = Intdc_pre + kp_v_dc*ydc_pre;
        Vq_pre = (Xg*Id_pre+Rg*Iq-Ug*sin(delta_pre)+Id_pre.*Lg.*Int_pre)./(1-Id_pre*Lg*kp_pll);
        omega_pre=kp_pll.*Vq_pre+Int_pre;
        t_prefault = [0;t_start];

        Iq_clear = Iq_pre(end);

        t_fault = [];
        delta_fault = [];
        Int_fault = [];
        ydc_fault = [];
        Intdc_fault = [];
        Vdc_fault = [];
        Id_fault = [];
        Vq_fault = [];
        Iq_fault = [];
        omega_fault=[];
      
        init = [prefault_SEP(1)-phasejump;prefault_SEP(2:5)];
        [t_postfault , x_all2] = ode78(@f_post,[t_start,t_end],init,odeset('RelTol',1e-6));
        delta_post = x_all2(:,1);
        Int_post = x_all2(:,2);
        ydc_post = x_all2(:,3);
        Intdc_post = x_all2(:,4);
        Iq_post = x_all2(:,5);
        Vdc_post = sqrt(ydc_post+Vdc_ref^2);
        Id_post = Intdc_post + kp_v_dc*ydc_post;
        Vq_post = (Xg*Id_post+Rg*Iq-Ug*sin(delta_post)+Id_post.*Lg.*Int_post)./(1-Id_post*Lg*kp_pll);
        omega_post=kp_pll.*Vq_post+Int_post;

         % fault-clear t0
        [t_fault_clear , x_all3] = ode78(@f_post,[t_start,t_start+t_0],init,odeset('RelTol',1e-10));
        delta_fault_clear = x_all3(:,1);
        Int_fault_clear = x_all3(:,2);
        ydc_fault_clear = x_all3(:,3);
        Intdc_fault_clear = x_all3(:,4);
        Iq_fault_clear = x_all3(:,5);
        Vdc_fault_clear = sqrt(ydc_fault_clear+Vdc_ref^2);
        Id_fault_clear = Intdc_fault_clear + kp_v_dc*ydc_fault_clear;
        Vq_fault_clear = (Xg*Id_fault_clear+Rg*Iq_fault_clear-Ug*sin(delta_fault_clear)+Id_fault_clear.*Lg.*Int_fault_clear)./(1-Id_fault_clear*Lg*kp_pll);
        omega_fault_clear=kp_pll.*Vq_fault_clear+Int_fault_clear;
            
        % reduce order
        if IsTVCslow == 1
            Iq = Iq_clear;  %for slow TVC
        else
            Iq = Iq00;
        end
        switch system
            case "DVCslow"
                delta_fault_RM2 =[];
                %postfault RM
                Intdc0 = Intdc_fault_clear(end);
                ydc0 = ydc_fault_clear(end);
                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[ydc0;Intdc0],odeset('RelTol',1e-10));
                Intdc_post_RM = x_all3(:,2);
                ydc_post_RM = x_all3(:,1);
                Id_post_RM = Intdc_post_RM + kp_v_dc*ydc_post_RM;
                Vdc_post_RM = sqrt(ydc_post_RM+Vdc_ref^2);
                delta_post_RM = asin((Xg*Id_post_RM+Rg*Iq)/Ug);

                % considering t0
                Intdc_post_RMc = [Intdc_fault_clear; Intdc_post_RM ];
                ydc_post_RMc = [ydc_fault_clear; ydc_post_RM];
                Id_post_RMc = Intdc_post_RMc + kp_v_dc*ydc_post_RMc;
                delta_post_RMc = asin((Xg*Id_post_RMc+Rg*Iq)/Ug);
                Vdc_post_RMc = sqrt(ydc_post_RMc+Vdc_ref^2);
                t_postfault_RMc = [t_fault_clear; t_postfault_RM];

             case "PLLslow"
                Id_fault_RM2 = [];
                %postfault RM
                Int0 = Int_fault_clear(end);
                delta0 = delta_fault_clear(end);

                Id_slow0 = (Pref+Iq.*Ug.*sin(delta0)-Iq.^2*Rg)./Ug./cos(delta0) - (Pref+Iq.*Ug.*sin(delta0)-Iq.^2*Rg).^2./Ug^3./cos(delta0).^3*Rg;
                ww = sqrt(Ug*cos(delta0)-Ug^2*cos(delta0)^2)/2;
                C1 = Id_fault_clear(end)-Id_slow0;
                C2 = (-C1*Ug*cos(delta0)/2+C_dc/2*zeta_dc*w_vdc*ydc_fault_clear(end))/ww;
                Id_star = (C1*Ug*cos(delta0)/2+C2*ww)/(Ug^2*cos(delta0)^2/4+ww^2);
                delta_star = 1/w_vdc*Id_star*(Xg*kp_pll)/(1-kp_pll*Lg*Id_slow0);
                Int_star = 1/w_vdc*Id_star*(Xg*ki_pll)/(1-kp_pll*Lg*Id_slow0);




                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[delta0+delta_star;Int0+Int_star],odeset('RelTol',1e-10));
                Int_post_RM = x_all3(:,2);
                delta_post_RM = x_all3(:,1);
                Id_post_RM = Pref./Ug./cos(delta_post_RM);

                % considering t0
                Int_post_RMc = [Int_fault_clear; Int_post_RM ];
                delta_post_RMc = [delta_fault_clear; delta_post_RM];
                Id_post_RMc = (Pref+Iq.*Ug.*sin(delta_post_RMc)-Iq.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq.*Ug.*sin(delta_post_RMc)-Iq.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
                t_postfault_RMc = [t_fault_clear; t_postfault_RM];


            case {"TVCfast", "PLLslow_TVC", "DVCslow_TVC"}
            clear x_all x_all2
            delta_fault_Rtvc = [];
            Int_fault_Rtvc = [];
            ydc_fault_Rtvc = [];
            Intdc_fault_Rtvc = [];   
            Vdc_fault_Rtvc = [];
            Id_fault_Rtvc = [];
            Iq_fault_Rtvc = [];
            Vq_fault_Rtvc = [];
            omega_fault_Rtvc=[];
            t_fault_Rtvc=[];
            

            init = [prefault_SEP(1)-phasejump;prefault_SEP(2:4)];
            [t_postfault_Rtvc , x_all2] = ode78(@f_TVCfast_post_time,[t_start,t_end],init,odeset('RelTol',1e-10));
            delta_post_Rtvc = x_all2(:,1);
            Int_post_Rtvc = x_all2(:,2);
            ydc_post_Rtvc = x_all2(:,3);
            Intdc_post_Rtvc = x_all2(:,4);
            Vdc_post_Rtvc = sqrt(ydc_post_Rtvc+Vdc_ref^2);
            Id_post_Rtvc = Intdc_post_Rtvc + kp_v_dc*ydc_post_Rtvc;
            Iq_post_Rtvc = (Ug*cos(delta_post_Rtvc)+Id_post_Rtvc*Rg-Vac_ref)/(Xg+1/m_v);
            Vq_post_Rtvc = (Xg*Id_post_Rtvc+Rg*Iq_post_Rtvc-Ug*sin(delta_post_Rtvc)+Id_post_Rtvc.*Lg.*Int_post_Rtvc)./(1-Id_post_Rtvc*Lg*kp_pll);
            omega_post_Rtvc=kp_pll.*Vq_post_Rtvc+Int_post_Rtvc;
        
            t_timedomain_Rtvc = [t_prefault;t_fault_Rtvc;t_postfault_Rtvc];
            delta_timedomain_Rtvc = [delta_pre; delta_fault_Rtvc; delta_post_Rtvc];
            Int_timedomain_Rtvc = [Int_pre; Int_fault_Rtvc; Int_post_Rtvc];%./2/pi+50;
            Vdc_timedomain_Rtvc = [Vdc_pre; Vdc_fault_Rtvc; Vdc_post_Rtvc];
            Id_timedomain_Rtvc = [Id_pre; Id_fault_Rtvc; Id_post_Rtvc];
            Iq_timedomain_Rtvc = [Iq_pre; Iq_fault_Rtvc; Iq_post_Rtvc];
            omega_timedomain_Rtvc = [omega_pre; omega_fault_Rtvc; omega_post_Rtvc];

             % fault-clear t0
            [t_fault_clear_Rtvc , x_all3] = ode78(@f_TVCfast_post_time,[t_start,t_start+t_0],init,odeset('RelTol',1e-10));
            delta_fault_clear_Rtvc = x_all3(:,1);
            Int_fault_clear_Rtvc = x_all3(:,2);
            ydc_fault_clear_Rtvc = x_all3(:,3);
            Intdc_fault_clear_Rtvc = x_all3(:,4);
            Vdc_fault_clear_Rtvc = sqrt(ydc_fault_clear_Rtvc+Vdc_ref^2);
            Id_fault_clear_Rtvc = Intdc_fault_clear_Rtvc + kp_v_dc*ydc_fault_clear_Rtvc;
            Iq_fault_clear_Rtvc = (Ug*cos(delta_fault_clear_Rtvc)+Id_fault_clear_Rtvc*Rg-Vac_ref)/(Xg+1/m_v);
            Vq_fault_clear_Rtvc = (Xg*Id_fault_clear_Rtvc+Rg*Iq_fault_clear_Rtvc-Ug*sin(delta_fault_clear_Rtvc)+Id_fault_clear_Rtvc.*Lg.*Int_fault_clear_Rtvc)./(1-Id_fault_clear_Rtvc*Lg*kp_pll);
            omega_fault_clear_Rtvc=kp_pll.*Vq_fault_clear_Rtvc+Int_fault_clear_Rtvc;

            if system == "DVCslow_TVC"   
                delta_fault_RM2 =[];
                Iq_fault_RM2 = [];
                %postfault RM
                Intdc0 = Intdc_fault_clear(end);
                ydc0 = ydc_fault_clear(end);
                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[ydc0;Intdc0],odeset('RelTol',1e-10));
                Intdc_post_RM = x_all3(:,2);
                ydc_post_RM = x_all3(:,1);
                Id_post_RM = Intdc_post_RM + kp_v_dc*ydc_post_RM;
                Vdc_post_RM = sqrt(ydc_post_RM+Vdc_ref^2);
                delta_post_RM = asin((Xg*Id_post_RM)/Ug);
                Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg - Vac_ref)*m_v/(Xg*m_v+1);
                delta_post_RM = asin((Id_post_RM*Xg+Iq_post_RM*Rg)/Ug);
                Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg - Vac_ref)*m_v/(Xg*m_v+1);


                % considering t0
                Intdc_post_RMc = [Intdc_fault_clear_Rtvc; Intdc_post_RM ];
                ydc_post_RMc = [ydc_fault_clear_Rtvc; ydc_post_RM];
                Id_post_RMc = Intdc_post_RMc + kp_v_dc*ydc_post_RMc;
                Vdc_post_RMc = sqrt(ydc_post_RMc+Vdc_ref^2);
                delta_post_RMc = asin((Xg*Id_post_RMc)/Ug);
                Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg - Vac_ref)*m_v/(Xg*m_v+1);
                delta_post_RMc = asin((Id_post_RMc*Xg+Iq_post_RMc*Rg)/Ug);
                Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg - Vac_ref)*m_v/(Xg*m_v+1);
                t_postfault_RMc = [t_fault_clear_Rtvc; t_postfault_RM];

            end
            if system == "PLLslow_TVC"
                Id_fault_RM2 = [];
                Iq_fault_RM2 = [];
                %postfault RM
                Int0 = Int_fault_clear(end);
                delta0 = delta_fault_clear(end);
                [t_postfault_RM, x_all3] = ode78(@f_reduce_post,[t_start+t_0,t_end],[delta0;Int0],odeset('RelTol',1e-10));
                Int_post_RM = x_all3(:,2);
                delta_post_RM = x_all3(:,1);
                Iq_post_RM = (Ug*cos(delta_post_RM) - Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RM = (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;
                Iq_post_RM = (Ug*cos(delta_post_RM) + Id_post_RM*Rg- Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RM = (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg)./Ug./cos(delta_post_RM) - (Pref+Iq_post_RM.*Ug.*sin(delta_post_RM)-Iq_post_RM.^2*Rg).^2./Ug^3./cos(delta_post_RM).^3*Rg;

                % considering t0
                Int_post_RMc = [Int_fault_clear_Rtvc; Int_post_RM ];
                delta_post_RMc = [delta_fault_clear_Rtvc; delta_post_RM];
                Iq_post_RMc = (Ug*cos(delta_post_RMc) - Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RMc = (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
                Iq_post_RMc = (Ug*cos(delta_post_RMc) + Id_post_RMc*Rg- Vac_ref)*m_v/(1+Xg*m_v);
                Id_post_RMc = (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg)./Ug./cos(delta_post_RMc) - (Pref+Iq_post_RMc.*Ug.*sin(delta_post_RMc)-Iq_post_RMc.^2*Rg).^2./Ug^3./cos(delta_post_RMc).^3*Rg;
                t_postfault_RMc = [t_fault_clear_Rtvc; t_postfault_RM];

            end
        end
end
%% draw the results
t_full_timedomain = [t_prefault;t_fault;t_postfault];
delta_timedomain = [delta_pre; delta_fault; delta_post];
Int_timedomain = [Int_pre; Int_fault; Int_post];%./2/pi+50;
Vdc_timedomain = [Vdc_pre; Vdc_fault; Vdc_post];
Id_timedomain = [Id_pre; Id_fault; Id_post];
Iq_timedomain = [Iq_pre; Iq_fault; Iq_post];
omega_timedomain = [omega_pre; omega_fault; omega_post];%./2/pi+50;

%delta -2
f2 = figure(2);
set(gcf,'position',[200 558 1300 300]);
grid on;hold on;
plot(t_full_timedomain,(mod(delta_timedomain+pi,2*pi)-pi)*180/pi,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);
%omega -3
f3 = figure(3);
set(gcf,'position',[200 558 1300 300]);
grid on;hold on;
plot(t_full_timedomain,omega_timedomain./2/pi+50,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);
% Int_delta -4
f4 = figure(4);
set(gcf,'position',[200 558 1300 300]);
grid on;hold on;
plot(t_full_timedomain,Int_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);
% Id -5
f5 = figure(5);
set(gcf,'position',[200 558 1300 300]);
grid on;hold on;
plot(t_full_timedomain,Id_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);
% Vdc -6
f6 = figure(6);
set(gcf,'position',[200 558 1300 300]);
grid on; hold on;
plot(t_full_timedomain,Vdc_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);
% Iq -7
f7 = figure(7);
set(gcf,'position',[200 558 1300 300]);
grid on; hold on;
plot(t_full_timedomain,Iq_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;     
yl=ylim;
ymin=yl(1,1);
ymax=yl(1,2);
trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
ylim([ymin,ymax]);
set(gca, 'FontSize', 14);


switch system
    case {"DVCslow", "DVCslow_TVC"} 
        t_RM2_timedomain = t_fault;
        %delta_RM2_timedomain = delta_fault_RM2;
        
        t_RM_timedomain = t_postfault_RMc;
        delta_RM_timedomain = delta_post_RMc;
        Vdc_RM_timedomain = Vdc_post_RMc;
        Id_RM_timedomain = Id_post_RMc;
    
        figure(f2);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        hold on;
        plot(t_RM_timedomain,delta_RM_timedomain*180/pi,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
        xlim([xmin,xmax]);

        figure(f3);
        hold on;
        plot(t_full_timedomain,(Xg*Id_timedomain-Ug.*sin(delta_timedomain))*kp_pll,'LineStyle','--','linewidth',2,'color','m');    hold on;
        diddt = kp_v_dc*2/(C_dc)*(Pref- (Id_post.*Ug.*cos(delta_post) + Id_post.^2*Rg))+ki_v_dc*ydc_post;
        plot(t_postfault,-1*(Xg*Id_post-Ug.*sin(delta_post))*kp_pll+Xg*diddt./sqrt(Ug^2-Xg*Id_post),'LineStyle','-','linewidth',2,'color','red');    hold on;

        
        figure(f5);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        hold on;
        plot(t_RM_timedomain,Id_RM_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
        xlim([xmin,xmax]);
        
        figure(f6);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        hold on;
        plot(t_RM_timedomain,Vdc_RM_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
        xlim([xmin,xmax]);       
        
        figure(f1);
        plot(Id_fault,ydc_fault,'k-','linewidth',1.5);
        plot(Id_post(1),ydc_post(1),'k.','MarkerSize',15);
        plot(Id_post,ydc_post,'k-','linewidth',1.5);
        
        plot(Id_post_RMc,ydc_post_RMc,'linewidth',1.5,'color',[0.8500 0.3250 0.0980]);

        if system == "DVCslow_TVC"
             Iq_RM2_timedomain = Iq_fault_RM2;
             Iq_RM_timedomain = Iq_post_RM;
             figure(f7);
             xl=xlim;
             xmin=xl(1,1);
             xmax=xl(1,2);
             hold on;
             plot(t_postfault_RM,Iq_RM_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
             plot(t_fault,Iq_RM2_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
             xlim([xmin,xmax]);
        end

    case {"PLLslow", "PLLslow_TVC"}
        
        t_RM2_timedomain = t_fault;
        
        t_RM_timedomain = t_postfault_RMc;
        delta_RM_timedomain = delta_post_RMc;
        Int_RM_timedomain = Int_post_RMc;
        Id_RM_timedomain = Id_post_RMc;
        
        figure(f2);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        hold on;
        plot(t_RM_timedomain,delta_RM_timedomain*180/pi,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
        xlim([xmin,xmax]); 

        figure(f3);
        hold on;
        plot(t_RM_timedomain,((Xg*Id_RM_timedomain-Ug.*sin(delta_RM_timedomain))*kp_pll+Int_RM_timedomain)./(1-kp_pll*Lg*Id_RM_timedomain)./2/pi+50,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;

        figure(f4);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        hold on;
        plot(t_RM_timedomain,Int_RM_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
        xlim([xmin,xmax]); 
        
   
        
        figure(f5);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        hold on;
        plot(t_RM_timedomain,Id_RM_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;

        
        xlim([xmin,xmax]); 
        
        figure(f1);
        plot(delta_fault,Int_fault,'k-','linewidth',1.5);
        plot(delta_post(1),Int_post(1),'k.','MarkerSize',15);
        plot(delta_post,Int_post,'k-','linewidth',1.5);
        
        plot(delta_post_RM,Int_post_RM,'linewidth',1.5,'color',[0.8500 0.3250 0.0980]);

        if system == "PLLslow_TVC"
             Iq_RM2_timedomain = Iq_fault_RM2;
             Iq_RM_timedomain = Iq_post_RM;
             figure(f7);
             xl=xlim;
             xmin=xl(1,1);
             xmax=xl(1,2);
             hold on;
             plot(t_postfault_RM,Iq_RM_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
             plot(t_fault,Iq_RM2_timedomain,'LineStyle','-','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
             xlim([xmin,xmax]);
        end

end

if (system == "TVCfast") || (system == "DVCslow_TVC" )  ||(system == "PLLslow_TVC"  )
        figure(f2);
        hold on;
        plot(t_timedomain_Rtvc,delta_timedomain_Rtvc*180/pi,'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
        
        figure(f3);
        hold on;
        plot(t_timedomain_Rtvc,omega_timedomain_Rtvc,'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);    
        
        figure(f4);
        hold on;
        plot(t_timedomain_Rtvc,Int_timedomain_Rtvc,'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);
        
        figure(f5);
        hold on;
        plot(t_timedomain_Rtvc,Id_timedomain_Rtvc,'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
        
        figure(f6);
        hold on;
        plot(t_timedomain_Rtvc,Vdc_timedomain_Rtvc,'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
       
    
        figure(f7);
        hold on;
        plot(t_timedomain_Rtvc,Iq_timedomain_Rtvc,'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;     
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);
end


save('initial.mat','Iq_clear','Iq_pre');% for slow TVC

%% function
function yes = isnewxep(ep_set,xep,torr)
    if isempty(ep_set)
        yes = 1;
        return;
    end
    minerr = inf;
    for m = 1 : length(ep_set)
        err = abs(xep - ep_set(m).xep);
        err = min(err, abs(2*pi-err));
        err = max(err);
        if minerr > err
            minerr = err;
        end
    end
    if(minerr>torr)
        yes = 1;
    else
        yes = 0;
    end
end



function dfdt = f(x)
global system;
switch system 
    case "DVCslow"
        dfdt = f_DVCslow(x);
    case "PLLslow"
        dfdt = f_PLLslow(x);
    case "overdamp"
        dfdt = f_overdamp(x);
    case "DVCslow_TVC"
        dfdt = f_DVCslowTVCfast(x);
    case "PLLslow_TVC"
        dfdt = f_PLLslowTVCfast(x);
    case "PLLfast_3D"  %test
        dfdt = f_PLLfast_3D(x);
end
end

function out = maxabs(in)

    out = abs(in);
    
    while length(out) > 1
        out = max(out);
    end

end

function disp_v(msg,v)
    disp([msg '=']);
    disp(v);
end

function dfdt = f_backward(t,x)
    dfdt = -f(x);
end

function dfdt = f_forward(t,x)
    dfdt = f(x);
end

function dfdt = f_reduce_fault(t,x)
global system;
switch system 
    case "DVCslow"
        dfdt = f_DVCslow_fault(x);
    case "DVCslow_TVC"
        dfdt = f_DVCslowTVCfast_fault(x);
    case "PLLslow_TVC"
        dfdt = f_PLLslowTVCfast_fault(x);
end
end
function dfdt = f_reduce_post(t,x)
global system;
switch system 
    case "DVCslow"
        dfdt = f_DVCslow_post(x);
    case "PLLslow"
        dfdt = f_PLLslow_post(x);
    case "DVCslow_TVC"
        dfdt = f_DVCslowTVCfast_post(x);
    case "PLLslow_TVC"
        dfdt = f_PLLslowTVCfast_post(x);
end
end
function dfdt = f_TVCfast_post_time(t,x)
    dfdt = f_TVCfast_post(x);
end
function dfdt = f_TVCfast_fault_time(t,x)
    dfdt = f_TVCfast_fault(x);
end
function dfdt = f_fault(t,x)
      dfdt = f_full_fault(x);
end
function dfdt = f_post(t,x)
      dfdt = f_full(x);
end
function dfdt = f_prefault(t,x)
      dfdt = f_full_pre(x);
end