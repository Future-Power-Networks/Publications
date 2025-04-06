% This function plots the phase protrait of theotical model, draw DOA, calculate SR 
% 

% Author(s): Yifan Zhang

%% parameter
    clc;
    I1d = 0.2;%0.8;
    I1q = 0;
    I2d = 0;
    I2q = 0;
    Ug = 1;%1

    k2p = 2.5*2*pi;

    Ws = 50*2*pi;
    
    R_ratio = 0.1;
    %post-fault
    X1 = 0.15;
    R1 = R_ratio *X1;
    X2 = 0.5;
    R2 = 0.005;
    Xg = 0.6;
    Rg = 0.04;

    %pre-fault
    X10 = X1;
    R10= R1;
    X20 = X2;
    R20= R2;
    Xg0 = Xg/2;
    Rg0 = Rg/2;

    Xeq = (X1*X2+X2*Xg+X1*Xg)/(X2+Xg);
    k1p = 2.5*2*pi/Xeq;

    %fault-on
    position = 0.8; %fault to inf bus
    Rf = 0.001;
    X1f = X1;
    R1f = R1;
    X2f = X2;
    R2f = R2;
    Im_temp = Rf*(Xg*position*1j+Rg*position)/(Rf + Xg*position*1j+Rg*position)+Xg*(1-position)*1j+Rg*(1-position);
    Imgf = Im_temp *(Xg*1j+Rg)/(Xg*1j+Rg+Im_temp); %(Rf//Xg*location+Xg*(1-location))//Xg
    Xg_f = imag(Imgf);
    Rg_f = real(Imgf);

%for GSP
    V2ref = 1;
    V1ref = 1;
    Kq = 4;

%for GFM
    Kgfm = 2.5*2*pi;
    Pref = 0.8;
    Vgfm = 1;

    Pref1 = I1d;
    Pref2 = Pref;
    Kgfm1 = k1p*Xeq;
    Kgfm2 = Kgfm;
    Vgfm1 = V1ref;
    Vgfm2 = Vgfm;


   
%for GFL with active power control
    Pgfl = 1;
    Kp_p = 0;
    Ki_p = 1*2*pi;

%fault set
    Iop_tmp = 1/((2-position)*(position)/2*(Xg*1j + Rg)+ Rf); %Xg*(2-location)//Xg*(location)+Rg
    Ug_fault = Iop_tmp*Rf+Iop_tmp/2*position*(1-position)*(Xg*1j + Rg);
    Ug_fault_angle = angle(Ug_fault);
    Ug_fault = abs(Ug_fault);
    
    t_c = 0.15; %0
    t_c2= 0.31; %2
    t_c3= 0.45; %4
    model = "rough";


    w_limit = 0.2*100*pi;

    figure;
    set(gcf,'position',[60 200 560 420]);

    hold on;
    grid on;

    CCR=[];
    X22222=[];
    Kqqqq=[];

%% function transient(N_machine)
global N; %#ok<*GVMIS>
global system; 
N=2;

system = "GFM";
state = 3;
switch state
    case 1
      system = "GFM2";
      Kgfm1 = 2.5*2*pi;
      color = [0 0 0];
    case 2
      system = "GFM2";
      Kgfm1 = 10*2*pi;
      color = [0 0 0];
    case 3
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 0;
      color = [0.94 0.8 0.1250];
    case 4
      system = "GFMGFLV2"; 
      k1p = 10*2*pi/Xeq;
      Kq = 0;
      color = [0.94 0.8 0.1250];
    case 5
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 1;
      color = [0.94 0.8 0.1250];
    case 6
      system = "GFMGFLV2"; 
      k1p = 10*2*pi/Xeq;
      Kq = 1;
      color = [0.94 0.8 0.1250];
    case 7
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 2;
      color = [0.9290 0.6940 0.1250];
    case 8
      system = "GFMGFLV2";
      k1p = 10*2*pi/Xeq;
      Kq = 2;
      color = [0.9290 0.6940 0.1250];
    case 9
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      color = [0.8500 0.3250 0.0980];
    case 10
      system = "GFMGFLV2";
      k1p = 10*2*pi/Xeq;
      Kq = 4;
      color = [0.8500 0.3250 0.0980];
    case 11
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 12;
      color = [0.6350 0.0780 0.1840];
    case 12
      system = "GFMGFLV2";
      k1p = 10*2*pi/Xeq;
      Kq = 12;
      color = [0.6350 0.0780 0.1840];

    case 13
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      Xg = 0.05;
      X2 = 0.95;
      color = [1 1 0.1250];
    case 14
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      Xg = 0.3;
      X2 = 0.8;
      color = [0.94 0.8 0.1250];
    case 15
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      Xg = 0.4;
      X2 = 0.7;
      color = [0.9290 0.6940 0.1250];
    case 16
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      Xg = 0.5;
      X2 = 0.6;
      color = [0.8500 0.3250 0.0980];
    case 17
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      Xg = 0.7;
      X2 = 0.4;
      color = [0.6350 0.0780 0.1840];
    case 18
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      Xg = 0.8;
      X2 = 0.3;
      color = [0.6350 0.0780 0.1840];
    case 19
      system = "GFMGFLV2";
      k1p = 2.5*2*pi/Xeq;
      Kq = 4;
      Xg = 0.95;
      X2 = 0.05;
      color = [0.6350 0.0780 0.1840];
      

end
for Kq = 0:0.5:4
    
torralence = 1e-2; 

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
%% calculate EPs
m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'display','off','MaxIterations',10000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);
    xep = mod(xep,2*pi);
    
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
            SEP = xep;
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
            prefault_sep = xep;
            end
        end
    end
end

for m = 1:length(ep_set_pre)
    disp_v('Index_pre',m);
    disp_v('Equilibrium_pre',ep_set_pre(m).xep);
    disp_v('Eigenvalue_pre', ep_set_pre(m).Lambda);
    disp_v('Eigenvector_pre',ep_set_pre(m).V);
end


%% 
%
if N == 2
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

    end

    color_code = {'blue','magenta','red','black'};
    if system == "GFMGFLP2"
         axis([0,4,-2*pi,2*pi])
    else
        axis([-1*pi,3/2*pi,-1*pi,3/2*pi]);
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
    

%     scatter(prefault_sep(1),prefault_sep(2),40,color_code{ep_set_pre(1).flag+1},'filled','LineWidth', 1.5);
    

    for m = 1 : length(ep_set_ext)
        xep = ep_set_ext(m).xep;
        flag= ep_set_ext(m).flag;
        scatter(xep(1),xep(2),60,color,'LineWidth', 1.5);
       
        if flag == 1
            v = ep_set_ext(m).v;
            perturb = 1e-3;
            switch system
                case "GFM"
                [~ , x_p] = ode45(@f_backward,[0,100],xep+v*perturb);
                [~ , x_n] = ode45(@f_backward,[0,100],xep-v*perturb);
                case "GFM2"
                [~ , x_p] = ode45(@f_backward,[0,100],xep+v*perturb);
                [~ , x_n] = ode45(@f_backward,[0,100],xep-v*perturb);
                case "STATCOM"
                [~ , x_p] = ode78(@f_backward,[0,30],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,30],xep-v*perturb,odeset('RelTol',1e-5));
                case "GFMGFL"
                [~ , x_p] = ode78(@f_backward,[0,60],xep+v*perturb,odeset('RelTol',1e-5));
                [~ , x_n] = ode78(@f_backward,[0,60],xep-v*perturb,odeset('RelTol',1e-5));
                case "GFMGFLV"
                [~ , x_p] = ode45(@f_backward,[0,30],xep+v*perturb);
                [~ , x_n] = ode45(@f_backward,[0,30],xep-v*perturb);
                case "GFMGFLV2"
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
            plot(x_all(:,1),x_all(:,2),'LineStyle','-','linewidth',2,'color',color);%scatter(x_all(:,1),x_all(:,2),'.');
            Stability_Radius = min([Stability_Radius,ad,bd]);
        end
        
    end
    %draw the stability circle
    %viscircles(SEP',Stability_Radius,'Color','b');%'Color','#A2142F')
    Stability_Radius
    CCR = [CCR Stability_Radius];
    X22222 = [X22222 X2];
    Kqqqq=[Kqqqq Kq];


%     %draw nullclines
%     x1=-2*pi:0.01*pi:2*pi;
%     x2=-2*pi:0.01*pi:2*pi;
%     [y1,y2]=meshgrid(x1,x2);
%     zz = zeros(length(x2),length(x1));
%     dzz = zeros(length(x2),length(x1));
%     nullcline1 = zeros(length(x2),length(x1));
%     nullcline2 = zeros(length(x2),length(x1));
%     for a = 1: length(x1)
%         for b = 1: length(x2)
%             switch system 
%             case "GFL"
%               vector = f2_GFL([y1(b,a);y2(b,a)]);
%             case "GFMGFL"
%               vector = f2_GFMGFL([y1(b,a);y2(b,a)]);
%             case "STATCOM"
%               vector = f2_STATCOM([y1(b,a);y2(b,a)]);
%             case "GFMGFLV"
%               vector = f2_GFMGFLV([y1(b,a);y2(b,a)]);
%             case "GFMGFLV2"
%               vector = f2_GFMGFLV2([y1(b,a);y2(b,a)]);
%             case "GFM"
%               vector = f2_GFM([y1(b,a);y2(b,a)]);
%             end
%             nullcline1(b,a) = vector(1);
%             nullcline2(b,a) = vector(2);
%         end
%     end
%     contour(y1,y2,nullcline1,[0 0],'g-','linewidth',1.5)
%     contour(y1,y2,nullcline2,[0 0],'y-','linewidth',1.5)

%     if system =="GFMGFLP2"
%     
%     x1= 0:0.1:5;
%     x2= -2*pi:0.15*pi:2*pi;
%     [y1,y2]=meshgrid(x1,x2);
%     yy1 = zeros(length(x2),length(x1));
%     yy2 = zeros(length(x2),length(x1));
%     for a = 1: length(x1)
%         for b = 1: length(x2)
%             yy = f2_GFMGFLP2([y1(b,a) y2(b,a)]);
%             yy1(b,a) = yy(1);
%             yy2(b,a) = yy(2);
%         end
%     end
%     quiver(y1,y2,yy1,yy2,'color',[0.6 0.6 0.6]);
%     xlabel('I_d');
%     ylabel('\delta_2/rad');
% 
%     else
% 
%     x1=-2*pi:0.15*pi:2*pi;
%     x2=-2*pi:0.15*pi:2*pi;
%     [y1,y2]=meshgrid(x1,x2);
%     yy1 = zeros(length(x1),length(x2));
%     yy2 = zeros(length(x1),length(x2));
%     for a = 1: length(x1)
%         for b = 1: length(x2)
%         switch system
%             case "GFM"
%             yy = f_GFM([y1(b,a) y2(b,a)]);
%             case "GFL"
%             yy = f2_GFL([y1(b,a) y2(b,a)]); 
%             case "STATCOM"
%             yy = f2_STATCOM([y1(b,a) y2(b,a)]);
%             case "GFMGFL"
%             yy = f2_GFMGFL([y1(b,a) y2(b,a)]);
%             case "GFMGFLV"
%             yy = f2_GFMGFLV([y1(b,a) y2(b,a)]);
%             case "GFMGFLP"
%             yy = f2_GFMGFLP([y1(b,a) y2(b,a)]);
%             case "General"
%             yy = f2_Generalized([y1(b,a) y2(b,a)]); 
%         end
%             yy1(b,a) = yy(1);
%             yy2(b,a) = yy(2);
%         end
%     end
%     %quiver(y1,y2,yy1,yy2,'color',[0.7 0.7 0.7]);
%     %xlabel('\delta_1/rad');
%     %ylabel('\delta_2/rad');
%     end



% %draw quiver
%     x1=-2*pi:0.15*pi:2*pi;
%     x2=-2*pi:0.15*pi:2*pi;
%     [y1,y2]=meshgrid(x1,x2);
%     yy1 = zeros(length(x1),length(x2));
%     yy2 = zeros(length(x1),length(x2));
%     for a = 1: length(x1)
%         for b = 1: length(x2)
%         switch system
%             case "GFM"
%             yy = f2_GFM([y1(b,a) y2(b,a)]);
%             case "GFL"
%             yy = f2_GFL([y1(b,a) y2(b,a)]); 
%             case "STATCOM"
%             yy = f2_STATCOM([y1(b,a) y2(b,a)]);
%             case "GFMGFL"
%             yy = f2_GFMGFL([y1(b,a) y2(b,a)]);
%             case "GFMGFLV"
%             yy = f2_GFMGFLV([y1(b,a) y2(b,a)]);
%             case "GFMGFLP"
%             yy = f2_GFMGFLP([y1(b,a) y2(b,a)]);
%             case "General"
%             yy = f2_Generalized([y1(b,a) y2(b,a)]); 
%         end
%             yy1(b,a) = yy(1);
%             yy2(b,a) = yy(2);
%         end
%     end
%     quiver(y1,y2,yy1,yy2,'color',[0.7 0.7 0.7]);


    %% draw fualt trajectory
%     [t1 , x_all] = ode78(@f_fault,[0,t_c],(prefault_sep-Ug_fault_angle),odeset('RelTol',1e-5));
%     x_all=x_all+Ug_fault_angle;
%     [tt , x_all2] = ode78(@f_forward,[0,2],[x_all(end,1),x_all(end,2)],odeset('RelTol',1e-5));
%     plot(x_all2(:,1),x_all2(:,2),'k','linewidth',1.5);
%     plot(x_all(end,1),x_all(end,2),'k.','MarkerSize',15);
%     plot(x_all2(1,1),x_all2(1,2),'k.','MarkerSize',15);
%     [t1f , x_allf] = ode78(@f_fault,[0,t_c2],prefault_sep-Ug_fault_angle,odeset('RelTol',1e-5));
%     x_allf=x_allf+Ug_fault_angle;
%     %plot(x_allf(:,1),x_allf(:,2),'k','linewidth',1.5)
%     [tt , x_all2] = ode78(@f_forward,[t1f(end),20],x_allf(end,:),odeset('RelTol',1e-5));
%     [t1f2 , x_allf2] = ode78(@f_fault,[0,t_c3],prefault_sep-Ug_fault_angle,odeset('RelTol',1e-5));
%     x_allf2=x_allf2+Ug_fault_angle;
%     plot(x_allf2(:,1),x_allf2(:,2),'k','linewidth',1.5);
%     [tt3 , x_all3] = ode78(@f_forward,[t1f2(end),20],x_allf2(end,:),odeset('RelTol',1e-5));
%     plot(x_all2(:,1),x_all2(:,2),'k','linewidth',1.5);
%     plot(x_allf(end,1),x_allf(end,2),'k.','MarkerSize',15);
%     plot(x_allf2(end,1),x_allf2(end,2),'k.','MarkerSize',15);
%     plot(x_all2(1,1),x_all2(1,2),'k.','MarkerSize',15);
%     plot(x_all3(:,1),x_all3(:,2),'k','linewidth',1.5);






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
end

function dfdt = f(x)
    global N; %#ok<GVMIS>
    global system;
    switch N
        case 2
            switch system 
                case "GFM"
                    dfdt = f2_GFM(x);
                case "GFM2"
                    dfdt = f2_GFM2(x);
                case "GFL"
                    dfdt = f2_GFL(x);
                case "STATCOM"
                    dfdt = f2_STATCOM(x);
                case "GFMGFL"
                    dfdt = f2_GFMGFL(x);
                case "GFMGFLV"
                    dfdt = f2_GFMGFLV(x);
                case "GFMGFLV2"
                    dfdt = f2_GFMGFLV2(x);
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
        case "GFMGFLV2"
          dfdt = f2_GFMGFLV_fault2(x);
        case "GFM"
          dfdt = f2_GFM_fault(x);
        case "GFM2"
          dfdt = f2_GFM2_fault(x);
    end
end
function dfdt = f_prefault(t,x)
global system;
    switch system 
        case "GFL"
          dfdt = f2_GFL_prefault(x);
        case "STATCOM"
          dfdt = f2_STATCOM_prefault(x);
        case "GFMGFL"
          dfdt = f2_GFMGFL_prefault(x);
        case "GFMGFLV"
          dfdt = f2_GFMGFLV_prefault(x);
        case "GFMGFLV2"
          dfdt = f2_GFMGFLV_prefault2(x);
        case "GFM"
          dfdt = f2_GFM_prefault(x);
        case "GFM2"
          dfdt = f2_GFM2_prefault(x);
    end
end
