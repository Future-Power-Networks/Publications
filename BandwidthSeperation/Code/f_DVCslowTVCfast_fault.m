function dfdt = f_DVCslowTVCfast_fault(x)

    Ug =evalin('base','Ug_fault');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');
    Id_limit = evalin('base','Id_limit');
    Lf = evalin('base','Lf');
    Vac_ref = evalin('base','Vac_ref');
    m_v = evalin('base','m_v');

    Lg = Xg/Ws;
    Lsum = Lg +Lf/Ws;

    y = x(1); %vdc^2-vdc_ref^2
    Int_id = x(2);

    fault_type = evalin('base','fault_type');
    R_model = evalin('base','R_model');
    
    if fault_type == "line_cut"
        Xg = evalin('base','Xg_f');
        Rg = evalin('base','Rg_f');
        Lg = Xg/Ws;
    end

    if Int_id >= Id_limit
        Int_id = Id_limit;
    end
    Id = Int_id + kp_v_dc*y;
    if Id >= Id_limit
        Id = Id_limit;
    end

    sindelta0 = (Id*Xg)/Ug;
    if sindelta0 >= 1
        sindelta0=1;
    end
    if sindelta0 <= -1
        sindelta0=-1;
    end
    cosdelta0 = sqrt(1-sindelta0^2);
    Iq0 = (Ug*cosdelta0 + Id*Rg - Vac_ref)*m_v/(Xg*m_v+1);
    sindelta = (Id*Xg+Iq0*Rg)/Ug;
    if sindelta >= 1
        sindelta=1;
    end
    if sindelta <= -1
        sindelta=-1;
    end
    cosdelta = sqrt(1-sindelta^2);
    Iq = (Ug*cosdelta + Id*Rg - Vac_ref)*m_v/(Xg*m_v+1);
    P = Id*Ug*cosdelta - Iq*Ug*sindelta + (Id^2+Iq^2)*Rg;  


    switch R_model
        case "normal"
            dfdt(1) = 2/(C_dc)*(Pref-P);% dy
            dfdt(2) = ki_v_dc*y; %d int_id
        case "precise"
            dfdt(1) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum)); %dy
            dfdt(2) = ki_v_dc*y; %dInt_id
    end

    dfdt = dfdt.';

    end