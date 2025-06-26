function dfdt = f_PLL_int(t,x)
        
    kp_pll = evalin('base','kp_pll');
    ki_pll = evalin('base','ki_pll');

    Iq = evalin('base','Iq00');  %for constant Iq
    Ug =evalin('base','Ug');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');
    Id_limit = evalin('base','Id_limit');
    Lf = evalin('base','Lf');

    Lg = Xg/Ws;
    Lsum = Lg +Lf/Ws;

    y = x(1); %vdc^2-vdc_ref^2
    Int_id = x(2);

    y_int = x(3);
    x_int = x(4);

    delta_sep = evalin('base','delta_pre');
    delta_sep = delta_sep(1);



    R_model = evalin('base','R_model');
    
    if Int_id >= Id_limit
        Int_id = Id_limit;
    end
    Id = Int_id + kp_v_dc*y;
    if Id >= Id_limit
        Id = Id_limit;
    end

    sindelta = (Id*Xg+Iq*Rg)/Ug;
    if sindelta >= 1
         sindelta = 1;
    elseif sindelta <=-1
         sindelta = -1;
    end
    cosdelta = sqrt(1-sindelta^2);

    delta = asin(sindelta);

    P = Id*Ug*cosdelta - Iq*Ug*sindelta + (Id^2+Iq^2)*Rg;

    da = ki_pll/kp_pll;

    aa = Xg^2*Id/(Ug^2-Xg^2*Id^2) / kp_pll;



    switch R_model
        case "rough"
            did = (ki_v_dc*y + kp_v_dc*2/(C_dc)*(Pref-P))/(1+kp_v_dc*2/C_dc*aa);
            dfdt(1) = 2/(C_dc)*(Pref-P-aa*did);% dy
            dfdt(2) = ki_v_dc*y; %d int_id
            
            dfdt(3) = -1*da*y_int-da^2*(delta-delta_sep);% dyint
            dfdt(4) = da/cosdelta/Ug*Xg*(ki_v_dc*y + kp_v_dc*2/(C_dc)*(Pref-P)) - da*x_int; %d xint
        case "normal"
            dfdt(1) = 2/(C_dc)*(Pref-P);% dy
            dfdt(2) = ki_v_dc*y; %d int_id
            dfdt(3) = -1*da*y_int-da^2*(delta-delta_sep);% dyint
            dfdt(4) = da/cosdelta/Ug*Xg*(ki_v_dc*y + kp_v_dc*2/(C_dc)*(Pref-P)) - da*x_int; %d xint
        case "precise"
            dfdt(1) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum)); %dy
            dfdt(2) = ki_v_dc*y; %dInt_id
            dfdt(3) = -1*da*y_int-da^2*(delta-delta_sep);% dyint
            dfdt(4) = da/cosdelta/Ug*Xg*(ki_v_dc*y + kp_v_dc*2/(C_dc)*(Pref-P)) - da*x_int; %d xint
    end

    dfdt = dfdt.';

end