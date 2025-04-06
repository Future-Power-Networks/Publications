% This function plots the phase protrait of Simulink and compares with therotical
% results

% Author(s): Yifan Zhang
%%
DeltaGFL1=log_Delta1.signals(1).values;
DeltaGFL2=log_Delta2.signals(1).values;


T_deta=Ts*10;
t_start= t0;
t_end1 = t_start+tc1;
t_end2 = 3;


DeltaGFL1_draw1=DeltaGFL1(t_start/T_deta:t_end1/T_deta+1);
DeltaGFL1_draw2=DeltaGFL1(t_end1/T_deta+1:t_end2/T_deta+1);
DeltaGFL2_draw1=DeltaGFL2(t_start/T_deta:t_end1/T_deta+1);
DeltaGFL2_draw2=DeltaGFL2(t_end1/T_deta+1:t_end2/T_deta+1);

for n = 1:length(DeltaGFL1_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) > 2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) + 2*pi;
    elseif (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL2_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL2_draw1(n)-DeltaGFL2_draw1(n+1)) > 2*pi*4/5
        DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) = DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) + 2*pi;
    elseif (DeltaGFL2_draw1(n)-DeltaGFL2_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) = DeltaGFL2_draw1(n+1:length(DeltaGFM_draw1)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL1_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) > 2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) + 2*pi;
    elseif (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL2_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL2_draw2(n)-DeltaGFL2_draw2(n+1)) > 2*pi*4/5
        DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) = DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) + 2*pi;
    elseif (DeltaGFL2_draw2(n)-DeltaGFL2_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) = DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) - 2*pi;
    end
end

%plot(DeltaGFL_draw1,DeltaGFM_draw1,'b','LineWidth',1.5);
hold on
plot(DeltaGFL2_draw2,DeltaGFL1_draw2,'b','LineWidth',1.5);
%save a.mat Theta_draw1 W_draw1;

DeltaGFL1=log_Delta1_2.signals(1).values;
DeltaGFL2=log_Delta2_2.signals(1).values;


T_deta=Ts*10;
t_start= t0;
t_end1 = t_start+tc2;
t_end2 = 3;


DeltaGFL1_draw1=DeltaGFL1(t_start/T_deta:t_end1/T_deta+1);
DeltaGFL1_draw2=DeltaGFL1(t_end1/T_deta+1:t_end2/T_deta+1);
DeltaGFL2_draw1=DeltaGFL2(t_start/T_deta:t_end1/T_deta+1);
DeltaGFL2_draw2=DeltaGFL2(t_end1/T_deta+1:t_end2/T_deta+1);
DeltaGFL2_draw2 =DeltaGFL2_draw2;



for n = 1:length(DeltaGFL1_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) > 2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) + 2*pi;
    elseif (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL2_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL2_draw1(n)-DeltaGFL2_draw1(n+1)) > 2*pi*4/5
        DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) = DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) + 2*pi;
    elseif (DeltaGFL2_draw1(n)-DeltaGFL2_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) = DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL1_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) > 2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) + 2*pi;
    elseif (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL2_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL2_draw2(n)-DeltaGFL2_draw2(n+1)) > 2*pi*4/5
        DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) = DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) + 2*pi;
    elseif (DeltaGFL2_draw2(n)-DeltaGFL2_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) = DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) - 2*pi;
    end
end
%plot(DeltaGFL1_draw1,DeltaGFL2_draw1,'r','LineWidth',1.5,'DisplayName','fault-on trajectories')
hold on
plot(DeltaGFL2_draw2,DeltaGFL1_draw2,'b','LineWidth',1.5,'DisplayName','post-fault trajectories')




DeltaGFL1=log_Delta1_3.signals(1).values;
DeltaGFL2=log_Delta2_3.signals(1).values;
T_deta=Ts*10;
t_start= t0;
t_end1 = t_start+tc3;
t_end2 = 3;
DeltaGFL1_draw1=DeltaGFL1(t_start/T_deta:t_end1/T_deta+1);
DeltaGFL1_draw2=DeltaGFL1(t_end1/T_deta+1:t_end2/T_deta+1);
DeltaGFL2_draw1=DeltaGFL2(t_start/T_deta:t_end1/T_deta+1);
DeltaGFL2_draw2=DeltaGFL2(t_end1/T_deta+1:t_end2/T_deta+1);
DeltaGFL2_draw2 =DeltaGFL2_draw2;



for n = 1:length(DeltaGFL1_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) > 2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) + 2*pi;
    elseif (DeltaGFL1_draw1(n)-DeltaGFL1_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) = DeltaGFL1_draw1(n+1:length(DeltaGFL1_draw1)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL2_draw1)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL2_draw1(n)-DeltaGFL2_draw1(n+1)) > 2*pi*4/5
        DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) = DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) + 2*pi;
    elseif (DeltaGFL2_draw1(n)-DeltaGFL2_draw1(n+1)) <  -2*pi*4/5
        DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) = DeltaGFL2_draw1(n+1:length(DeltaGFL2_draw1)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL1_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) > 2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) + 2*pi;
    elseif (DeltaGFL1_draw2(n)-DeltaGFL1_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) = DeltaGFL1_draw2(n+1:length(DeltaGFL1_draw2)) - 2*pi;
    end
end
for n = 1:length(DeltaGFL2_draw2)-1  % Check the "continuous" property of phase angle
    if (DeltaGFL2_draw2(n)-DeltaGFL2_draw2(n+1)) > 2*pi*4/5
        DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) = DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) + 2*pi;
    elseif (DeltaGFL2_draw2(n)-DeltaGFL2_draw2(n+1)) <  -2*pi*4/5
        DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) = DeltaGFL2_draw2(n+1:length(DeltaGFL2_draw2)) - 2*pi;
    end
end
plot(DeltaGFL2_draw1,DeltaGFL1_draw1,'r','LineWidth',1.5,'DisplayName','fault-on trajectories')
hold on
plot(DeltaGFL2_draw2,DeltaGFL1_draw2,'b','LineWidth',1.5,'DisplayName','post-fault trajectories')




