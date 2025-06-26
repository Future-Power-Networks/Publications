function dfdt = f_GFL(x)

    Id = evalin('base','Id');
    Iq = evalin('base','Iq');
    Ug =evalin('base','Ug');
        
    kp = evalin('base','kp_pll');
    ki = evalin('base','ki_pll');

    Lg = evalin('base','Lg');
    Xg = evalin('base','Xg');


    delta = x(1);
    Int = x(2);

    Vq = (Xg*Id-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp);

    dfdt(1) = kp*Vq+Int;
    dfdt(2) = ki*Vq; %dint


    dfdt = dfdt.';

    end