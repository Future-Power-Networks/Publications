%%

DeltaGFL=ScopeData.signals(6).values; %degree
OmegaGFL=(ScopeData.signals(5).values-1)*Wbase; %rad
IdGFL=ScopeData.signals(3).values(:,1); %Id
IqGFL=ScopeData.signals(3).values(:,2); %Id
VdcGFL=ScopeData.signals(4).values; %rad
ydcGFL = VdcGFL.^2-Vdc_ref.^2;
IntGFL = ScopeData1.signals(1).values;
IntdcGFL = ScopeData1.signals(2).values;


t_GFL = ScopeData.time;


T_deta=Ts*10;
t_start2 = t_sim_start+t_start;
t_end1 = t_start2+t_c;
t_end2 = t_sim_start+t_end;

global system
switch system
    case {"DVCslow","DVCslow_TVC"}
        X_draw1=IdGFL(t_start2/T_deta+1:t_end1/T_deta+1);
        X_draw2=IdGFL(t_end1/T_deta+1:t_end2/T_deta+1);
        Y_draw1=ydcGFL(t_start2/T_deta+1:t_end1/T_deta+1);
        Y_draw2=ydcGFL(t_end1/T_deta+1:t_end2/T_deta+1);
    case {"PLLslow","PLLslow_TVC"}
        X_draw1=DeltaGFL(t_start2/T_deta+1:t_end1/T_deta+1)/180*pi;
        X_draw2=DeltaGFL(t_end1/T_deta+1:t_end2/T_deta+1)/180*pi;
        Y_draw1=IntGFL(t_start2/T_deta+1:t_end1/T_deta+1);
        Y_draw2=IntGFL(t_end1/T_deta+1:t_end2/T_deta+1);
end

try
figure(f1);
plot(X_draw1,Y_draw1,'r-','LineWidth',1.5,'DisplayName','fault-on trajectories')
hold on
plot(X_draw2,Y_draw2,'r-','LineWidth',1.5,'DisplayName','post-fault trajectories')
catch
end

delta_simulation = DeltaGFL(t_sim_start/T_deta+1:t_end2/T_deta+1)/180*pi;
omega_simulation = OmegaGFL(t_sim_start/T_deta+1:t_end2/T_deta+1);
Int_simulation = IntGFL(t_sim_start/T_deta+1:t_end2/T_deta+1);
Intdc_simulation = IntdcGFL(t_sim_start/T_deta+1:t_end2/T_deta+1);
Id_simulation = IdGFL(t_sim_start/T_deta+1:t_end2/T_deta+1);
Vdc_simulation = VdcGFL(t_sim_start/T_deta+1:t_end2/T_deta+1);
t_simulation = t_GFL(t_sim_start/T_deta+1:t_end2/T_deta+1) -t_sim_start;
Iq_simulation = IqGFL(t_sim_start/T_deta+1:t_end2/T_deta+1);


%%
figure(f2);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
plot(t_simulation,delta_simulation*180/pi,'LineStyle','-','linewidth',2,'color',[0 0.4470 0.7410]);    hold on;
xlim([xmin,xmax]);



figure(f3);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
plot(t_simulation,omega_simulation./2/pi+50,'LineStyle','-','linewidth',2,'color',[0 0.4470 0.7410]);    hold on;
xlim([xmin,xmax]);


figure(f4);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
plot(t_simulation,Int_simulation,'LineStyle','-','linewidth',2,'color',[0 0.4470 0.7410]);    hold on;
xlim([xmin,xmax]);

figure(f5);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
plot(t_simulation,Id_simulation,'LineStyle','-','linewidth',2,'color',[0 0.4470 0.7410]);    hold on;
xlim([xmin,xmax]);

figure(f6);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on; 
plot(t_simulation,Vdc_simulation,'LineStyle','-','linewidth',2,'color',[0 0.4470 0.7410]);    hold on;
xlim([xmin,xmax]);


figure(f7);
xl=xlim;
xmin=xl(1,1);
xmax=xl(1,2);
hold on;
plot(t_simulation,Iq_simulation,'LineStyle','-','linewidth',2,'color',[0 0.4470 0.7410]);    hold on;
xlim([xmin,xmax]);

savefig(f1,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSimMC1'));
savefig(f2,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSimMC2'));
savefig(f3,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSimMC3'));
savefig(f4,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSimMC4'));
savefig(f5,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSimMC5'));
savefig(f6,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSimMC6'));
savefig(f7,strcat('C:\Users\yz7521\OneDrive - Imperial College London\Desktop\DSP\FigureSimMC7'));
%%
experiment = 1;
save(strcat('C:\Users\yz7521\OneDrive - Imperial College London\MCIB\test',num2str(experiment),'.mat'),'t_end');

