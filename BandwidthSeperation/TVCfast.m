    
%% this function for evaluate the TVC loop which is configured as vary fast 
    close all

    global model; 
    model = "precise2";%"precise" "normal"
    global system;
    system  = "TVCfast";
    global R_model; 
    R_model = "normal2";%"normal2";%"num1" "normal" "precise"

%%
    Iq = 0;
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
    Vq_fault = (Xg*Id_fault+Rg*Iq-Ug_fault*sin(delta_fault)+Id_fault.*Lg.*Int_fault)./(1-Id_fault*Lg*kp_pll);
    omega_fault=kp_pll.*Vq_fault+Int_fault;

    dIq_precise_fault = zeros(size(x_all,1),1);
    ddelta_precise_fault = zeros(size(x_all,1),1);
    dy_precise_fault = zeros(size(x_all,1),1);
    dIq_normal_fault = zeros(size(x_all,1),1);
    ddelta_normal_fault = zeros(size(x_all,1),1);
    dy_normal_fault = zeros(size(x_all,1),1);
    for i=1:size(x_all,1)
        model = "precise2";
        dfdttemp = f_full_fault(x_all(i,:));
        dIq_precise_fault(i) = dfdttemp(5);
        ddelta_precise_fault(i) = dfdttemp(1);
        dy_precise_fault(i) = dfdttemp(3);
        model = "normal2";
        dfdttemp = f_full_fault(x_all(i,:));
        dIq_normal_fault(i) = dfdttemp(5);
        ddelta_normal_fault(i) = dfdttemp(1);
        dy_normal_fault(i) = dfdttemp(3);
        model = "precise2";
    end

  

    [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,:),odeset('RelTol',1e-10));
    delta_post = x_all2(:,1);
    Int_post = x_all2(:,2);
    ydc_post = x_all2(:,3);
    Intdc_post = x_all2(:,4);
    Iq_post = x_all2(:,5);
    Vdc_post = sqrt(ydc_post+Vdc_ref^2);
    Id_post = Intdc_post + kp_v_dc*ydc_post;
    Vq_post = (Xg*Id_post+Rg*Iq-Ug*sin(delta_post)+Id_post.*Lg.*Int_post)./(1-Id_post*Lg*kp_pll);
    omega_post=kp_pll.*Vq_post+Int_post;

    dIq_precise_post = zeros(size(x_all2,1),1);
    ddelta_precise_post = zeros(size(x_all2,1),1);
    dy_precise_post = zeros(size(x_all2,1),1);
    dIq_normal_post = zeros(size(x_all2,1),1);
    ddelta_normal_post = zeros(size(x_all2,1),1);
    dy_normal_post = zeros(size(x_all,1),1);
    for i=1:size(x_all2,1)
        model = "precise2";
        dfdttemp = f_full(x_all2(i,:));
        dIq_precise_post(i) = dfdttemp(5);
        ddelta_precise_post(i) = dfdttemp(1);
        dy_precise_post(i) = dfdttemp(3);
        model = "normal2";
        dfdttemp = f_full(x_all2(i,:));
        dIq_normal_post(i) = dfdttemp(5);
        ddelta_normal_post(i) = dfdttemp(1);
        dy_normal_post(i) = dfdttemp(3);
        model = "precise2";
    end
   

    t_full_timedomain = [t_prefault;t_fault;t_postfault];
    delta_timedomain = [delta_pre; delta_fault; delta_post];
    Int_timedomain = [Int_pre; Int_fault; Int_post];%./2/pi+50;
    Vdc_timedomain = [Vdc_pre; Vdc_fault; Vdc_post];
    Id_timedomain = [Id_pre; Id_fault; Id_post];
    Iq_timedomain = [Iq_pre; Iq_fault; Iq_post];






    %%
    
        f2 = figure(2);
        set(gcf,'position',[200 558 1300 300]);
        grid on;hold on;
        plot(t_full_timedomain,delta_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
       
       
        %plot(t_RM1_timedomain,delta_RM1_timedomain,'LineStyle','-','linewidth',2,'color',[255/255 0/255 95/255]);    hold on;
        
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);
        
        f3 = figure(3);
        set(gcf,'position',[200 558 1300 300]);
        grid on;hold on;
        plot(t_full_timedomain,Int_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);
        
        f4 = figure(4);
        set(gcf,'position',[200 558 1300 300]);
        grid on;hold on;
        plot(t_full_timedomain,Id_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
        
        
        %plot(t_RM1_timedomain,Id_RM1_timedomain,'LineStyle','-','linewidth',2,'color',[255/255 0/255 95/255]);    hold on;
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);
        
        f5 = figure(5);
        set(gcf,'position',[200 558 1300 300]);
        grid on; hold on;
        plot(t_full_timedomain,Vdc_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
        %plot(t_RM1_timedomain,Vdc_RM1_timedomain,'LineStyle','-','linewidth',2,'color',[255/255 0/255 95/255]);    hold on;
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);

        f6 = figure(6);
        set(gcf,'position',[200 558 1300 300]);
        grid on; hold on;
        plot(t_full_timedomain,Iq_timedomain,'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;     
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);
    


%         figure(7)
%         set(gcf,'position',[200 558 1300 300]);
%         grid on; hold on;
%         plot([t_fault;t_postfault],[dIq_precise_fault;dIq_precise_post],'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;  
%         plot([t_fault;t_postfault],[dIq_normal_fault;dIq_normal_post],'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);
%         plot([t_fault;t_postfault],[dIq_precise_fault;dIq_precise_post]-[dIq_normal_fault;dIq_normal_post],'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
% 
%         figure(8)
%         set(gcf,'position',[200 558 1300 300]);
%         grid on; hold on;
%         plot([t_fault;t_postfault],[ddelta_precise_fault;ddelta_precise_post],'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;  
%         plot([t_fault;t_postfault],[ddelta_normal_fault;ddelta_normal_post],'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);
%         plot([t_fault;t_postfault],[ddelta_precise_fault;ddelta_precise_post]-[ddelta_normal_fault;ddelta_normal_post],'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
% 
%         figure(9)
%         set(gcf,'position',[200 558 1300 300]);
%         grid on; hold on;
%         plot([t_fault;t_postfault],[dy_precise_fault;dy_precise_post],'LineStyle','-','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;  
%         plot([t_fault;t_postfault],[dy_normal_fault;dy_normal_post],'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);
%         plot([t_fault;t_postfault],[dy_precise_fault;dy_precise_post]-[dy_normal_fault;dy_normal_post],'LineStyle','-','linewidth',2,'color',[0.4940 0.1840 0.5560]);
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
%%
    model = "normal2";
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
    Vq_fault = (Xg*Id_fault+Rg*Iq-Ug_fault*sin(delta_fault)+Id_fault.*Lg.*Int_fault)./(1-Id_fault*Lg*kp_pll);
    omega_fault=kp_pll.*Vq_fault+Int_fault;
  

    [t_postfault , x_all2] = ode78(@f_post,[t_fault(end),t_end],x_all(end,:),odeset('RelTol',1e-10));
    delta_post = x_all2(:,1);
    Int_post = x_all2(:,2);
    ydc_post = x_all2(:,3);
    Intdc_post = x_all2(:,4);
    Iq_post = x_all2(:,5);
    Vdc_post = sqrt(ydc_post+Vdc_ref^2);
    Id_post = Intdc_post + kp_v_dc*ydc_post;
    Vq_post = (Xg*Id_post+Rg*Iq-Ug*sin(delta_post)+Id_post.*Lg.*Int_post)./(1-Id_post*Lg*kp_pll);
    omega_post=kp_pll.*Vq_post+Int_post;

    t_full_timedomain = [t_prefault;t_fault;t_postfault];
    delta_timedomain = [delta_pre; delta_fault; delta_post];
    Int_timedomain = [Int_pre; Int_fault; Int_post];%./2/pi+50;
    Vdc_timedomain = [Vdc_pre; Vdc_fault; Vdc_post];
    Id_timedomain = [Id_pre; Id_fault; Id_post];
    Iq_timedomain = [Iq_pre; Iq_fault; Iq_post];

    figure(f2);
    set(gcf,'position',[200 558 1300 300]);
    grid on;hold on;
    plot(t_full_timedomain,delta_timedomain,'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    
   
    %plot(t_RM1_timedomain,delta_RM1_timedomain,'LineStyle','-','linewidth',2,'color',[255/255 0/255 95/255]);    hold on;
    
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    
    figure(f3);
    set(gcf,'position',[200 558 1300 300]);
    grid on;hold on;
    plot(t_full_timedomain,Int_timedomain,'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    
    figure(f4);
    set(gcf,'position',[200 558 1300 300]);
    grid on;hold on;
    plot(t_full_timedomain,Id_timedomain,'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
   
    
    %plot(t_RM1_timedomain,Id_RM1_timedomain,'LineStyle','-','linewidth',2,'color',[255/255 0/255 95/255]);    hold on;
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    
    figure(f5);
    set(gcf,'position',[200 558 1300 300]);
    grid on; hold on;
    plot(t_full_timedomain,Vdc_timedomain,'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
   
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);

    figure(f6);
    set(gcf,'position',[200 558 1300 300]);
    grid on; hold on;
    plot(t_full_timedomain,Iq_timedomain,'LineStyle','-','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;     
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);



    %% reduce model
    clear x_all x_all2

    [t_fault_R , x_all] = ode78(@f_reduce_fault,[t_start,t_start+t_c],prefault_SEP(1:4),odeset('RelTol',1e-10));
    delta_fault_R = x_all(:,1);
    Int_fault_R = x_all(:,2);
    ydc_fault_R = x_all(:,3);
    Intdc_fault_R = x_all(:,4);   
    Vdc_fault_R = sqrt(ydc_fault_R+Vdc_ref^2);
    Id_fault_R = Intdc_fault_R + kp_v_dc*ydc_fault_R;
    Iq_fault_R = (Ug_fault*cos(delta_fault_R)+Id_fault_R*Rg-Vac_ref)/(Xg+1/m_v);
    Vq_fault_R = (Xg*Id_fault_R+Rg*Iq_fault_R-Ug_fault*sin(delta_fault_R)+Id_fault_R.*Lg.*Int_fault_R)./(1-Id_fault_R*Lg*kp_pll);
    omega_fault_R=kp_pll.*Vq_fault_R+Int_fault_R;
    
  

    [t_postfault_R , x_all2] = ode78(@f_reduce_post,[t_fault_R(end),t_end],x_all(end,:),odeset('RelTol',1e-10));
    delta_post_R = x_all2(:,1);
    Int_post_R = x_all2(:,2);
    ydc_post_R = x_all2(:,3);
    Intdc_post_R = x_all2(:,4);
    Vdc_post_R = sqrt(ydc_post_R+Vdc_ref^2);
    Id_post_R = Intdc_post_R + kp_v_dc*ydc_post_R;
    Iq_post_R = (Ug*cos(delta_post_R)+Id_post_R*Rg-Vac_ref)/(Xg+1/m_v);
    Vq_post_R = (Xg*Id_post_R+Rg*Iq_post_R-Ug*sin(delta_post_R)+Id_post_R.*Lg.*Int_post_R)./(1-Id_post_R*Lg*kp_pll);
    omega_post_R=kp_pll.*Vq_post_R+Int_post_R;

    t_timedomain_R = [t_prefault;t_fault_R;t_postfault_R];
    delta_timedomain_R = [delta_pre; delta_fault_R; delta_post_R];
    Int_timedomain_R = [Int_pre; Int_fault_R; Int_post_R];%./2/pi+50;
    Vdc_timedomain_R = [Vdc_pre; Vdc_fault_R; Vdc_post_R];
    Id_timedomain_R = [Id_pre; Id_fault_R; Id_post_R];
    Iq_timedomain_R = [Iq_pre; Iq_fault_R; Iq_post_R];

    figure(f2);
    set(gcf,'position',[200 558 1300 300]);
    grid on;hold on;
    plot(t_timedomain_R,delta_timedomain_R,'LineStyle','-','linewidth',2,'color','green');    hold on;
    
   
    %plot(t_RM1_timedomain,delta_RM1_timedomain,'LineStyle','-','linewidth',2,'color',[255/255 0/255 95/255]);    hold on;
    
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    
    figure(f3);
    set(gcf,'position',[200 558 1300 300]);
    grid on;hold on;
    plot(t_timedomain_R,Int_timedomain_R,'LineStyle','-','linewidth',2,'color','green');    hold on;
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    
    figure(f4);
    set(gcf,'position',[200 558 1300 300]);
    grid on;hold on;
    plot(t_timedomain_R,Id_timedomain_R,'LineStyle','-','linewidth',2,'color','green');    hold on;
   
    
    %plot(t_RM1_timedomain,Id_RM1_timedomain,'LineStyle','-','linewidth',2,'color',[255/255 0/255 95/255]);    hold on;
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    
    figure(f5);
    set(gcf,'position',[200 558 1300 300]);
    grid on; hold on;
    plot(t_timedomain_R,Vdc_timedomain_R,'LineStyle','-','linewidth',2,'color','green');    hold on;
   
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);

    figure(f6);
    set(gcf,'position',[200 558 1300 300]);
    grid on; hold on;
    plot(t_timedomain_R,Iq_timedomain_R,'LineStyle','-','linewidth',2,'color','green');    hold on;     
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    trange=[t_start,t_start+t_c,t_start+t_c,t_start];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);

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
    case "TVCfast"
        dfdt = f_TVCfast_fault(x);
end
end
function dfdt = f_reduce_post(t,x)
global system;
switch system 
    case "DVCslow"
        dfdt = f_DVCslow_post(x);
    case "PLLslow"
        dfdt = f_PLLslow_post(x);
    case "TVCfast"
        dfdt = f_TVCfast_post(x);
end
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