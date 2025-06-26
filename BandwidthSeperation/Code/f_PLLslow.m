function dfdt = f_PLLslow(x)


    Iq = evalin('base','Iq');  %for constant Iq

    Ug =evalin('base','Ug');
        
    kp_pll = evalin('base','kp_pll');
    ki_pll = evalin('base','ki_pll');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    Lg = Xg/Ws;
    Pref = evalin('base','Pref');

    delta = x(1);
    Int = x(2);
    if delta >= (pi/2 - 1e-3)
        delta = (pi/2 - 1e-3);
    elseif delta <= (-pi/2 + 1e-3)
        delta = (-pi/2 + 1e-3);
    end

    Id = (Pref+Iq*Ug*sin(delta)-Iq^2*Rg)/Ug/cos(delta) - (Pref+Iq*Ug*sin(delta)-Iq^2*Rg)^2/Ug^3/cos(delta)^3*Rg; %tylor expression


    Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp_pll);

    if delta >= (pi/2 - 1e-3)
        dfdt(1) = 0; %ddelta
        dfdt(2) = 0; %dint
    elseif delta <= (-pi/2 + 1e-3)
        dfdt(1) = 0; %ddelta
        dfdt(2) = 0; %dint
    else
        dfdt(1) = kp_pll*Vq+Int; %ddelta
        dfdt(2) = ki_pll*Vq; %dint
    end


    dfdt = dfdt.';

    end