function dfdt = f_overdamp(x)

    Iq = evalin('base','Iq');
    Ug =evalin('base','Ug');
        
    kp_v_dc = evalin('base','kp_v_dc');

    kp_pll = evalin('base','kp_pll');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    Lg = Xg/Ws;
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');

    R_model = evalin('base','R_model');

    Lf = evalin('base','Lf');
    Lf = Lf / Ws;
    Lsum = Lf + Lg;

    delta = x(1);
    Id = x(2);


    switch R_model
        case "normal"
            P = Id*Ug*cos(delta) - Iq*Ug*sin(delta) + (Id^2+Iq^2)*Rg;
            Vq = (Xg*Id+Rg*Iq-Ug*sin(delta))/(1-Id*Lg*kp_pll);

            dfdt(1) = kp_pll*Vq; %ddelta
            dfdt(2) = kp_v_dc*(Pref - P)*2/C_dc; %d id

        case "precise"
            P = Id*Ug*cos(delta) - Iq*Ug*sin(delta) + (Id^2+Iq^2)*Rg;
            Vq = (Xg*Id+Rg*Iq-Ug*sin(delta))/(1-Id*Lg*kp_pll);

            dfdt(1) = kp_pll*Vq; %ddelta
            dfdt(2) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref- P); %did
    end





    dfdt = dfdt.';

    end