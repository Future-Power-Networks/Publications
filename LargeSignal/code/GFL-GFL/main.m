% This function plots the phase protrait of theotical model, draw DOA, calculate SR 
% 

% Author(s): Yifan Zhang


%% parameter

    I1d = 0.8;%0.8;
    I1q = 0;
    I2d = 0.4;%0.8;
    I2q = 0;
    Ug = 1;%1

    k1p = 10*2*pi;%2.5*2*pi;
    k2p = 10*2*pi;

    Ws = 50*2*pi;
    


    %post-fault
    X1 = 0.2;%0.12;
    X2 = 0.2;%0.12;
    Xg = 0.7;%0.7;

    %pre-fault
    X10 = X1;
    X20 = X2;
    Xg0 = Xg/2;

    %fault-on
    position = 0.5; %fault to inf bus
    Rf = 0.02;
    X1f = X1;
    X2f = X2;
    Im_temp = Rf*(Xg*position*1j)/(Rf + Xg*position*1j)+Xg*(1-position)*1j;
    Imgf = Im_temp *Xg*1j/(Xg*1j+Im_temp); %(Rf//Xg*location+Xg*(1-location))//Xg
    Xg_f = imag(Imgf);
    Rg_f = real(Imgf);




%for gfl with voltage control
    V1ref = 1;
    Kq = 2;

%for GFM
    Kgfm = 1*2*pi;
    Vgfm = 1;%1.0309
    Pref = 0.6;
   
%for GFL with active power control
    Pgfl = 1;
    Kp_p = 0;
    Ki_p = 1*2*pi;

%fault set
    Iop_tmp = 1/((2-position)*(position)/2*Xg*1j + Rf); %Xg*(2-location)//Xg*(location)+Rg
    Ug_fault = Iop_tmp*Rf+Iop_tmp/2*position*(1-position)*Xg*1j;
    Ug_fault_angle = angle(Ug_fault);
    Ug_fault= abs(Ug_fault);
    
    t_c = 0.212;
    t_c2=0.212;
    model = "rough";%"resistance";%"rough"; %"resistance";  %"precise"; %"rough" %"resistance"


    w_limit = 0.2*100*pi;




   % model = "precise"; %"rough";  %"precise"; "rough_limit" "precise_limit"



%% function transient(N_machine)
global N; %#ok<*GVMIS>
global system; 
N=2;
system = "GFL";
%for Xg=0.7:0.02:0.7

torralence = 1e-2; 

try
    N;
catch
    N = 2;
end

try
    system;
catch
    system = "GFL";
end

if system == "GFMGFLP2"
x=(0:0.1:1)*5;
y=(0:0.1:1)*2*pi;
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
else
x=(0:0.1:1)*2*pi;
n = length(x);
x_set = zeros(N,n^N);
for m = 0:(n^N - 1)
    for l = 1:N
        index = floor(m/n^(l-1));
        index = mod(index,n)+1;
        x_set(l,m+1) = x(index);
    end
end
end
%% calculate EPs
m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',10000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);
    if system == "GFMGFLP2"
        xep(2) = mod(real(xep(2)),2*pi);
    else
        xep = mod(xep,2*pi);
    end
    
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
clear xep flag v V Lambda A sig m
ep_set_pre = [];
m = 1;
%prefault
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@(x)f_prefault(0,x),xep,options);

   
    if maxabs(ferr) < torralence
        if isnewxep(ep_set_pre,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set_pre(m).xep = real(xep); %#ok<*SAGROW> 
            ep_set_pre(m).A = A;
            ep_set_pre(m).Lambda = Lambda;
            ep_set_pre(m).V = V;   
            ep_set_pre(m).v = v;     % stable eigenvectors of unstable ep 
            ep_set_pre(m).flag = flag;
           
            m = m+1;
            if flag == 0
            jacob=A;
            xeps=xep;
            end
        end
    end
end

for m = 1:length(ep_set)
    disp_v('Index_pre',m);
    disp_v('Equilibrium_pre',ep_set_pre(m).xep);
    disp_v('Eigenvalue_pre', ep_set_pre(m).Lambda);
    disp_v('Eigenvector_pre',ep_set_pre(m).V);
end


%% 
%
if N == 2
if system == "GFMGFLP2"
    clear ep_set_ext;
     for n = 1:length(ep_set)
        mm = (n-1)*3;
        ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(mm+2)=ep_set(n);
        ep_set_ext(mm+3)=ep_set(n); 
        ep_set_ext(mm+2).xep(2) = ep_set(n).xep(2) - 2*pi;
        ep_set_ext(mm+3).xep(2) = ep_set(n).xep(2) + 2*pi;
     end
else
    clear ep_set_ext;
    for n = 1:length(ep_set)
        m = (n-1)*6;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
        ep_set_ext(m+4)=ep_set(n);
        ep_set_ext(m+5)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+6)=ep_set(n);
        ep_set_ext(m+7)=ep_set(n);
        ep_set_ext(m+8)=ep_set(n);
    
        ep_set_ext(m+2).xep(1) = ep_set(n).xep(1) - 2*pi;
    
        ep_set_ext(m+3).xep(2) = ep_set(n).xep(2) - 2*pi;
    
        ep_set_ext(m+4).xep(1) = ep_set(n).xep(1) - 2*pi;
        ep_set_ext(m+4).xep(2) = ep_set(n).xep(2) - 2*pi;
        ep_set_ext(m+5).xep(1) = ep_set(n).xep(1);
        ep_set_ext(m+5).xep(2) = ep_set(n).xep(2) - 4*pi;

        ep_set_ext(m+6).xep(1) = ep_set(n).xep(1) - 2*pi;
        ep_set_ext(m+6).xep(2) = ep_set(n).xep(2) - 4*pi;
%         ep_set_ext(m+6).xep(1) = ep_set(n).xep(1) - 6*pi;
%         ep_set_ext(m+6).xep(2) = ep_set(n).xep(2) - 2*pi;
%         ep_set_ext(m+7).xep(1) = ep_set(n).xep(1) - 4*pi;
%         ep_set_ext(m+7).xep(2) = ep_set(n).xep(2);
%         ep_set_ext(m+8).xep(1) = ep_set(n).xep(1) - 6*pi;
%         ep_set_ext(m+8).xep(2) = ep_set(n).xep(2);
%         ep_set_ext(m+9).xep(1) = ep_set(n).xep(1) - 8*pi;
%         ep_set_ext(m+9).xep(2) = ep_set(n).xep(2);
%         ep_set_ext(m+10).xep(1) = ep_set(n).xep(1) - 10*pi;
%         ep_set_ext(m+10).xep(2) = ep_set(n).xep(2);
%         ep_set_ext(m+11).xep(1) = ep_set(n).xep(1) - 12*pi;
%         ep_set_ext(m+11).xep(2) = ep_set(n).xep(2);
%         ep_set_ext(m+12).xep(1) = ep_set(n).xep(1) - 14*pi;
%         ep_set_ext(m+12).xep(2) = ep_set(n).xep(2);

    end
end    
    figure;
    %title((['X_g is ',num2str(Xg),' Xg']))
    hold on;
    grid on;
    color_code = {'blue','magenta','red','black'};
    if system == "GFMGFLP2"
         axis([0,4,-2*pi,2*pi])
    else
        axis([-pi,3/2*pi,-pi,3/2*pi]);
%         xticks(-3/2*pi:pi/8:3/2*pi);
%         yticks(-3/2*pi:pi/8:3/2*pi);
%         xticklabels({'$-1\frac{1}{2}\pi$', '', '$-1\frac{1}{4}\pi$', '', '$-1\pi$', '', '$-\frac{3}{4}\pi$','', '$-\frac{1}{2}\pi$','','$-\frac{1}{4}\pi$','',...
%              '0','', '$\frac{1}{4}\pi$', '', '$\frac{1}{2}\pi$', '', '$\frac{3}{4}\pi$', '', '$\pi$', '', ...
%              '$1\frac{1}{4}\pi$', '', '$1\frac{1}{2}\pi$', '', '$1\frac{3}{2}\pi$'});
%         yticklabels({'$-1\frac{1}{2}\pi$', '', '$-1\frac{1}{4}\pi$', '', '$-1\pi$', '', '$-\frac{3}{4}\pi$','', '$-\frac{1}{2}\pi$','','$-\frac{1}{4}\pi$','',...
%              '0','', '$\frac{1}{4}\pi$', '', '$\frac{1}{2}\pi$', '', '$\frac{3}{4}\pi$', '', '$\pi$', '', ...
%              '$1\frac{1}{4}\pi$', '', '$1\frac{1}{2}\pi$', '', '$1\frac{3}{2}\pi$'});
         xticks(-2*pi:pi/2:2*pi);
         yticks(-2*pi:pi/2:2*pi);
         xticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
         yticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
         set(gca, 'TickLabelInterpreter', 'latex');
         set(gca, 'FontSize', 14);
    end
    %model = "rough_limit";
    %SEP = ep_set(find([ep_set.flag]==0)).xep;
    DISTANCE = [];
    ad = 0;
    bd = 0;
    Stability_Radius = 2*pi;
    prefault_sep=ep_set_pre(1).xep;

    scatter(prefault_sep(1),prefault_sep(2),40,color_code{ep_set_pre(1).flag+1},'filled','LineWidth', 1.5);
    SEP = ep_set(1).xep;

    for m = 1 : length(ep_set_ext)
        xep = ep_set_ext(m).xep;
        flag= ep_set_ext(m).flag;
        scatter(xep(1),xep(2),60,color_code{flag+1},'LineWidth', 1.5);
       
        if flag == 1
            v = ep_set_ext(m).v;
            perturb = 1e-3;
            switch system
                case "GFM"
                [~ , x_p] = ode45(@f_backward,[0,1000],xep+v*perturb);
                [~ , x_n] = ode45(@f_backward,[0,1000],xep-v*perturb);
                case "STATCOM"
                [~ , x_p] = ode45(@f_backward,[0,30],xep+v*perturb);
                [~ , x_n] = ode45(@f_backward,[0,30],xep-v*perturb);
                case "GFMGFL"
                [~ , x_p] = ode45(@f_backward,[0,60],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode45(@f_backward,[0,60],xep-v*perturb,odeset('RelTol',1e-5));
                case "GFMGFLV"
                [~ , x_p] = ode45(@f_backward,[0,30],xep+v*perturb);
                [~ , x_n] = ode45(@f_backward,[0,30],xep-v*perturb);
                case "GFMGFLP"
                [~ , x_p] = ode45(@f_backward,[0,30],xep+v*perturb);
                [~ , x_n] = ode45(@f_backward,[0,30],xep-v*perturb);
                case "GFMGFLP2"
                [~ , x_p] = ode45(@f_backward,[0,200],xep+v*perturb*100,odeset('RelTol',1e-3));
                [~ , x_n] = ode45(@f_backward,[0,200],xep-v*perturb*100,odeset('RelTol',1e-3));
                case "GFL"
                [~ , x_p] = ode45(@f_backward,[0,60],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode45(@f_backward,[0,60],xep-v*perturb,odeset('RelTol',1e-5));
                case "General"
                [~ , x_p] = ode45(@f_backward,[0,50],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode45(@f_backward,[0,50],xep-v*perturb,odeset('RelTol',1e-5));
            end

            DISTANCE = x_p-ones(length(x_p),1)*SEP';
            ad = min(hypot(DISTANCE(:,1),DISTANCE(:,2)));
            DISTANCE = x_n-ones(length(x_n),1)*SEP';
            bd = min(hypot(DISTANCE(:,1),DISTANCE(:,2)));
            x_all = [flip(x_n,1);x_p];
            plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);%scatter(x_all(:,1),x_all(:,2),'.');
            Stability_Radius = min([Stability_Radius,ad,bd]);
        end
        
    end
    %draw the stability circle
    viscircles(SEP',Stability_Radius,'Color','b');%'Color','#A2142F')
    Stability_Radius

    if system =="GFMGFLP2"
    
    x1= 0:0.1:5;
    x2= -2*pi:0.15*pi:2*pi;
    [y1,y2]=meshgrid(x1,x2);
    yy1 = zeros(length(x2),length(x1));
    yy2 = zeros(length(x2),length(x1));
    for a = 1: length(x1)
        for b = 1: length(x2)
            yy = f2_GFMGFLP2([y1(b,a) y2(b,a)]);
            yy1(b,a) = yy(1);
            yy2(b,a) = yy(2);
        end
    end
    quiver(y1,y2,yy1,yy2,'color',[0.6 0.6 0.6]);
    xlabel('I_d');
    ylabel('\delta_2/rad');

    else

    x1=-2*pi:0.15*pi:2*pi;
    x2=-2*pi:0.15*pi:2*pi;
    [y1,y2]=meshgrid(x1,x2);
    yy1 = zeros(length(x1),length(x2));
    yy2 = zeros(length(x1),length(x2));
    for a = 1: length(x1)
        for b = 1: length(x2)
        switch system
            case "GFM"
            yy = f_GFM([y1(b,a) y2(b,a)]);
            case "GFL"
            yy = f2_GFL([y1(b,a) y2(b,a)]); 
            case "STATCOM"
            yy = f2_STATCOM([y1(b,a) y2(b,a)]);
            case "GFMGFL"
            yy = f2_GFMGFL([y1(b,a) y2(b,a)]);
            case "GFMGFLV"
            yy = f2_GFMGFLV([y1(b,a) y2(b,a)]);
            case "GFMGFLP"
            yy = f2_GFMGFLP([y1(b,a) y2(b,a)]);
            case "General"
            yy = f2_Generalized([y1(b,a) y2(b,a)]); 
        end
            yy1(b,a) = yy(1);
            yy2(b,a) = yy(2);
        end
    end
    %quiver(y1,y2,yy1,yy2,'color',[0.7 0.7 0.7]);
    %xlabel('\delta_1/rad');
    %ylabel('\delta_2/rad');
    end


    %% draw fualt trajectory
    [t1 , x_all] = ode78(@f_fault,[0,t_c],(prefault_sep-Ug_fault_angle),odeset('RelTol',1e-5));
    x_all=x_all+Ug_fault_angle;
    [tt , x_all2] = ode78(@f_forward,[0,1],[x_all(end,1),x_all(end,2)],odeset('RelTol',1e-5));
    plot(x_all2(:,1),x_all2(:,2),'k','linewidth',1.5);
    plot(x_all(end,1),x_all(end,2),'k.','MarkerSize',15);
    plot(x_all2(1,1),x_all2(1,2),'k.','MarkerSize',15);
    [t1 , x_all] = ode78(@f_fault,[0,t_c2],prefault_sep-Ug_fault_angle,odeset('RelTol',1e-5));
    x_all=x_all+Ug_fault_angle;
    plot(x_all(:,1),x_all(:,2),'k','linewidth',1.5)
    [tt , x_all2] = ode78(@f_forward,[t1(end),20],x_all(end,:),odeset('RelTol',1e-5));
    plot(x_all2(:,1),x_all2(:,2),'k','linewidth',1.5);
    plot(x_all(end,1),x_all(end,2),'k.','MarkerSize',15);
    plot(x_all2(1,1),x_all2(1,2),'k.','MarkerSize',15);

   [t1 , x_all] = ode78(@f_forward,[0,1],[-1*pi-4*pi,0],odeset('RelTol',1e-5));
   plot(x_all(:,1),x_all(:,2),'g-','MarkerSize',15);



% three apparatus
elseif N == 3

    figure;
    color_code = {'blue','magenta','green','red','cyan','blue','magenta','green','red','cyan'};
    for m = 1:length(ep_set)
        xep = ep_set(m).xep;
        flag = ep_set(m).flag;
        scatter3(xep(1),xep(2),xep(3),color_code{flag+1});
        hold on;
    end    

    ep_set_ext = ep_set;
    n = 1;
    for m = 1 : length(ep_set_ext)        
        flag = ep_set_ext(m).flag;
        if flag 
            xep = ep_set_ext(m).xep;
            v = ep_set_ext(m).v;
            perturb = 10e-2;

            if flag == 1
                for alpha = (0:0.01:1)*2*pi
                    vp = v(:,1)*sin(alpha) + v(:,2)*cos(alpha);
                    [~ , x_all] = ode45(@f_backward,[0,1000],xep+vp*perturb);
                    scatter3(x_all(:,1),x_all(:,2),x_all(:,3),'.',color_code{n});
                end
                n = n + 1;
            elseif flag == 2
                for beta = [-1,1]
                    vp = v*beta;
                    [~ , x_all] = ode45(@f_backward,[0,1000],xep+vp*perturb);
                    scatter3(x_all(:,1),x_all(:,2),x_all(:,3),'.',color_code{n});
                end
                n = n + 1;
            end
            
        end
    end

end
%end
%end

function dfdt = f(x)
    global N; %#ok<GVMIS>
    global system;
    switch N
        case 2
            switch system 
                case "GFM"
                    dfdt = f_GFM(x);
                case "GFL"
                    dfdt = f2_GFL(x);
                case "STATCOM"
                    dfdt = f2_STATCOM(x);
                case "GFMGFL"
                    dfdt = f2_GFMGFL(x);
                case "GFMGFLV"
                    dfdt = f2_GFMGFLV(x);
                case "GFMGFLP"
                    dfdt = f2_GFMGFLP(x);
                case "General"
                    dfdt = f2_Generalized(x);
                case "GFMGFLP2"
                    dfdt = f2_GFMGFLP2(x);
            end
        case 3
            dfdt = f3(x);
    end
end

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

function disp_v(msg,v)
    disp([msg '=']);
    disp(v);
end

function dfdt = f_forward(t,x)
    dfdt = f(x);
end

function dfdt = f_backward(t,x)
    dfdt = -f(x);
end

function dfdt = f_fault(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f2_GFL_fault(x);
        case "STATCOM"  
          dfdt = f2_STATCOM_fault(x); 
        case "GFMGFL"
          dfdt = f2_GFMGFL_fault(x);
        case "GFMGFLP2"
          dfdt = f2_GFMGFLP2_fault(x);
        case "GFMGFLV"
          dfdt = f2_GFMGFLV_fault(x);
    end
end
function dfdt = f_prefault(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f2_GFL_prefault(x);
    end
end
