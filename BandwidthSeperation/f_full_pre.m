function dfdt = f_full_pre(x)

    Iq0 = evalin('base','Iq');
    Ug =evalin('base','Ug');
        
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
    Vac_ref = evalin('base','Vac_ref');
    w_tvc = evalin('base','w_tvc');
    m_v = evalin('base','m_v');


    fault_type = evalin('base','fault_type');

    if fault_type == "line_cut"
        Xg0 = evalin('base','Xg0');
        Rg0 = evalin('base','Rg0');
        Xg = Xg0;
        Rg = Rg0;
        Lg = Xg/Ws;
    elseif fault_type == "power_jump"
        Pref = evalin('base','Pref0');
    end

    delta = x(1);
    Int = x(2);
    y = x(3); %vdc^2-Vdc_ref^2
    Vdc = sqrt(y+Vdc_ref^2);
    Int_id = x(4);

    Iq_f = x(5); %iq low-pass filter state
    Iq = Iq_f + Iq0;
    
    Id = Int_id + kp_v_dc*y;
    P = Id*Ug*cos(delta) - Iq*Ug*sin(delta) + (Id^2+Iq^2)*Rg;
    Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp_pll);
    Vd = Ug*cos(delta) - Iq*Xg + Id*Rg;


    dfdt(1) = kp_pll*Vq+Int; %ddelta
    dfdt(2) = ki_pll*Vq; %dint

    dfdt(3) = 2/C_dc*(Pref-P); %dy
    dfdt(4) = ki_v_dc*y; %dInt_id

    dfdt(5) = w_tvc*(-Iq_f + m_v*(Vd-Vac_ref) );  %d Iq_f

    dfdt = dfdt.';

    end