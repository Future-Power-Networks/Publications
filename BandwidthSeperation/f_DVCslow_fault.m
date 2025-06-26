function dfdt = f_DVCslow_fault(x)

    Iq = evalin('base','Iq');
    Ug =evalin('base','Ug_fault');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');
    Id_limit = evalin('base','Id_limit');

    y = x(1); %vdc^2-vdc_ref^2
    Int_id = x(2);
    if Int_id >= Id_limit
        Int_id = Id_limit;
    end
    Id = Int_id + kp_v_dc*y;
    if Id >= Id_limit
        Id = Id_limit;
    end

    sindelta = Id*Xg/Ug;
    cosdelta = sqrt(1-sindelta^2);
    if (Id<Ug/Xg) && (Id>-Ug/Xg)
        P = Id*Ug*cosdelta - Iq*Ug*sindelta + (Id^2+Iq^2)*Rg;
    else
        P=0;
    end

    R_model = evalin('base','R_model');

    switch R_model
        case "rough"
            dfdt(1) = 2/(C_dc)*(Pref-P);% dy
            dfdt(2) = ki_v_dc*y; %d int_id
        case "normal"
            dfdt(1) = 2/(C_dc)*(Pref-P);% dy
            dfdt(2) = ki_v_dc*y; %d int_id
        case "precise"
            dfdt(1) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum)); %dy
            dfdt(2) = ki_v_dc*y; %dInt_id
    end

    dfdt(1) = 2/(C_dc)*(Pref-P);% dy
    dfdt(2) = ki_v_dc*y; %d int_id

    dfdt = dfdt.';
    end