function dfdt = f_TVCfast_post(x)
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

    Id_limit = evalin('base','Id_limit');


    R_model = evalin('base','R_model_TVC');

    w_tvc = evalin('base','w_tvc');
    m_v = evalin('base','m_v');
    Vac_ref = evalin('base','Vac_ref');

    Lf = evalin('base','Lf');
    Lf = Lf / Ws;
    Lsum = Lf + Lg;
    
    y = x(3); %vdc^2-Vdc_ref^2


    Vdc = sqrt(y+Vdc_ref^2);
    Int_id = x(4);
    if Int_id >= Id_limit
        Int_id = Id_limit;
    end
    Id = Int_id + kp_v_dc*y;
    if Id >= Id_limit
        Id = Id_limit;
    end

    delta = x(1);
    Int = x(2);


    switch R_model
        case "normal"
            Iq = (Ug*cos(delta) + Id*Rg - Vac_ref)/(Xg+1/m_v);


            P = Id*Ug*cos(delta) - Iq*Ug*sin(delta) + (Id^2+Iq^2)*Rg;
            Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp_pll);

            dfdt(1) = kp_pll*Vq+Int; %ddelta
            dfdt(2) = ki_pll*Vq; %dint
      
            dfdt(3) = 2/C_dc*(Pref-P); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id

        case "normal2"

            Iq = (Ug*cos(delta) + Id*Rg - Vac_ref)/(Xg+1/m_v);


            P = Id*Ug*cos(delta) - Iq*Ug*sin(delta) + (Id^2+Iq^2)*Rg;
            Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp_pll);

            dfdt(1) = kp_pll*Vq+Int; %ddelta
            dfdt(2) = ki_pll*Vq; %dint
      
            dfdt(3) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum )); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id
            
        case "precise"

            Iq = (Ug*cos(delta) + Id*Rg - Vac_ref)/(Xg+1/m_v);


            P = Id*Ug*cos(delta) - Iq*Ug*sin(delta) + (Id^2+Iq^2)*Rg;
            Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp_pll);

            dfdt(1) = (kp_pll*Vq+Int)/(1+Lg*Ug*sin(delta)*kp_pll/(1-kp_pll*Lg*Id)); %ddelta
            dfdt(2) = ki_pll*Vq; %dint
      
            dfdt(3) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum)); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id

      
    end

    dfdt = dfdt.';

    end