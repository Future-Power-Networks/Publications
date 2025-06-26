function f = Fun_AEfslove_DVCslow_fault(x,ydc,Intdc)

    Iq = evalin('base','Iq');
    Ug =evalin('base','Ug_fault');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    kp_pll = evalin('base','kp_pll');
    ki_pll = evalin('base','ki_pll');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    Lg = Xg/Ws;
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');

    Vdc_ref= evalin('base','Vdc_ref');

    
    
    y = ydc; %vdc^2-Vdc_ref^2
    if y<=-Vdc_ref^2
        y=-Vdc_ref^2;
    end
    Vdc = sqrt(y+Vdc_ref^2);
    Int_id = Intdc;
    Id = Int_id + kp_v_dc*y;
    delta = x(1);
    Int = x(2);
    
    Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp_pll);


    f(1) = kp_pll*Vq+Int; %ddelta
    f(2) = ki_pll*Vq; %dint

    end