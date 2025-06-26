function dfdt = f_DVCslow_post(x)


    Iq = evalin('base','Iq');  %for constant Iq
    Ug =evalin('base','Ug');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');
    Id_limit = evalin('base','Id_limit');
    kp_pll = evalin('base','kp_pll');
    ki_pll = evalin('base','ki_pll');
    Lf = evalin('base','Lf');

    SEP = evalin('base','prefault_SEP');

    Lg = Xg/Ws;
    Lsum = Lg +Lf/Ws;

    y = x(1); %vdc^2-vdc_ref^2
    Int_id = x(2);

    R_model = evalin('base','R_model');
    
    if Int_id >= Id_limit
        Int_id = Id_limit;
    end
    Id = Int_id + kp_v_dc*y;
    if Id >= Id_limit
        Id = Id_limit;
    end

    sindelta = (Id*Xg+Iq*Rg)/Ug;
    cosdelta = sqrt(1-sindelta^2);
    if (Id<Ug/Xg) && (Id>-Ug/Xg)
        P = Id*Ug*cosdelta - Iq*Ug*sindelta + (Id^2+Iq^2)*Rg;
    else
        P=0;
    end

    ipsi = Lg*Id;


    switch R_model
        case "rough"
            dfdt(1) = 2/(C_dc)*(Pref-P);% dy
            dfdt(2) = ki_v_dc*y; %d int_id
        case "normal"
            dfdt(1) = 2/(C_dc)*(Pref-P);% dy
            dfdt(2) = ki_v_dc*y; %d int_id
        case "num"
            a = Xg^2*Id^2/(Ug^2-Xg^2*Id^2)*(1/kp_pll-Id*Lg);
            b = Xg*Id/sqrt(Ug^2-Xg^2*Id^2)*ki_pll/kp_pll^2*(asin(sindelta)-SEP(1));
            dfdt(1) = 2/(C_dc+2*kp_v_dc*Id*Lsum +2*a*kp_v_dc)*(Pref-(P + ki_v_dc*y*Id*Lsum + a*ki_v_dc*y-b)); %dy
            dfdt(2) = ki_v_dc*y; %dInt_id
        case "precise"
            dfdt(1) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum)); %dy
            dfdt(2) = ki_v_dc*y; %dInt_id
        case "num_kp"
            a = Xg^2*Id^2/(Ug^2-Xg^2*Id^2)*(1/kp_pll-Id*Lg);
            dfdt(1) = 2/(C_dc+2*kp_v_dc*Id*Lsum +2*a*kp_v_dc)*(Pref-(P + ki_v_dc*y*Id*Lsum + a*ki_v_dc*y)); %dy
            dfdt(2) = ki_v_dc*y; %dInt_id
    end

    dfdt = dfdt.';

    end