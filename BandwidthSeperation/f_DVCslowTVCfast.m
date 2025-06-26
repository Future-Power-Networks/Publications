function dfdt = f_DVCslowTVCfast(x)

    Ug =evalin('base','Ug');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');
    m_v = evalin('base','m_v');
    Vac_ref = evalin('base','Vac_ref');
    R_model = evalin('base','R_model');
    Lf = evalin('base','Lf');
    Lg = Xg/Ws;
    Lf = Lf / Ws;
    Lsum = Lf + Lg;
    Id_limit = evalin('base','Id_limit');
    

    Id = x(1);
    y = x(2); %vdc^2-vdc_ref^2

    if (Id<Ug/Xg) && (Id>-Ug/Xg)
        sindelta0 = (Id*Xg)/Ug;
        cosdelta0 = sqrt(1-sindelta0^2);
        Iq0 = (Ug*cosdelta0 + Id*Rg - Vac_ref)*m_v/(Xg*m_v+1);
        sindelta = (Id*Xg+Iq0*Rg)/Ug;
        cosdelta = sqrt(1-sindelta^2);
        Iq = (Ug*cosdelta + Id*Rg - Vac_ref)*m_v/(Xg*m_v+1);
        P = Id*Ug*cosdelta - Iq*Ug*sindelta + (Id^2+Iq^2)*Rg;
    else
        P=0;
    end

    if Id>Id_limit

        Id =Id_limit;
    end


    switch R_model
    case "normal"
        dfdt(1) = ki_v_dc*y+kp_v_dc*2/(C_dc)*(Pref-P);   %dId
        dfdt(2) = 2/(C_dc)*(Pref-P);  %dy
    case "precise"
        dfdt(1) = ki_v_dc*y+kp_v_dc*2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref- (P+ ki_v_dc*y*Id*Lsum));   %dId
        dfdt(2) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref- (P+ ki_v_dc*y*Id*Lsum));  %dy    
    end
    dfdt = dfdt.';

    end