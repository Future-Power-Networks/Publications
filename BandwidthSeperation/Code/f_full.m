function dfdt = f_full(x)

    Iq0 = evalin('base','Iq00');
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

    model = evalin('base','model');

    w_tvc = evalin('base','w_tvc');
    m_v = evalin('base','m_v');
    Vac_ref = evalin('base','Vac_ref');

    Lf = evalin('base','Lf');
    Lf = Lf / Ws;
    Lsum = Lf + Lg;

    
    
    y = x(3); %vdc^2-Vdc_ref^2

    Iq_f = x(5); %iq low-pass filter state
    Iq = Iq_f + Iq0;

    if y<=-Vdc_ref^2
        y=-Vdc_ref^2;
    end
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
    Vd = Ug*cos(delta) - Iq*Xg + Id*Rg;

    w_tvc_eq = w_tvc/( 1 + w_tvc*m_v*Iq*Lg^2*kp_pll/(1-Id*Lg*kp_pll) + w_tvc*m_v*Lg*kp_v_dc*2/C_dc*Iq*Lsum/(1+2/C_dc*kp_v_dc*Id*Lsum) );
    P = Id*Ug*cos(delta) - Iq*Ug*sin(delta) + (Id^2+Iq^2)*Rg;
    Vq = (Xg*Id+Rg*Iq-Ug*sin(delta)+Id*Lg*Int)/(1-Id*Lg*kp_pll);

    switch model
        case "rough"
            Vq = Xg*Id+Rg*Iq-Ug*sin(delta);
            diqdt = w_tvc*(-Iq_f + m_v*(Vd-Vac_ref));   %d Iq_f

            dfdt(1) = kp_pll*Vq+Int; %ddelta
            dfdt(2) = ki_pll*Vq; %dint
      
            dfdt(3) = 2/C_dc*(Pref-P); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id

            dfdt(5) = diqdt;  %d Iq_f

        case "normal"
            
            diqdt = w_tvc*(-Iq_f + m_v*(Vd-Vac_ref));   %d Iq_f

            dfdt(1) = kp_pll*Vq+Int; %ddelta
            dfdt(2) = ki_pll*Vq; %dint
      
            dfdt(3) = 2/C_dc*(Pref-P); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id

            dfdt(5) = diqdt;  %d Iq_f
        
        case "normal2"

            diqdt = w_tvc_eq*(-Iq_f + m_v*(Vd-Vac_ref) );   %d Iq_f

            dfdt(1) = kp_pll*Vq+Int; %ddelta
            dfdt(2) = ki_pll*Vq; %dint
      
            dfdt(3) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum)); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id

            dfdt(5) = diqdt;  %d Iq_f

        case "precise2"

            diqdt = w_tvc_eq*(-Iq_f + m_v*(Vd-Vac_ref) - m_v*Iq*Lg*(kp_pll*Vq+Int) + m_v*Lg*(ki_v_dc*y + kp_v_dc*2/C_dc*(Pref-P))/(1+2/C_dc*Id*kp_v_dc*Lsum) );   %d Iq_f

            dfdt(1) = kp_pll*Vq+Int + kp_pll*Lg*diqdt/(1-Id*Lg*kp_pll); %ddelta
            dfdt(2) = ki_pll*Vq + ki_pll*Lg*diqdt/(1-Id*Lg*kp_pll); %dint
      
            dfdt(3) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum + diqdt*Iq*Lsum )); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id

            dfdt(5) = diqdt;  %d Iq_f


        case "precise"

            diqdt = w_tvc*(-Iq_f + m_v*(Vd-Vac_ref) );   %d Iq_f

            dfdt(1) = kp_pll*Vq + Int + Lg*diqdt/(1-Id*Lg*kp_pll); %ddelta
            dfdt(2) = ki_pll*Vq; %dint
      
            dfdt(3) = 2/(C_dc+2*kp_v_dc*Id*Lsum)*(Pref-(P + ki_v_dc*y*Id*Lsum + diqdt*Iq*Lsum)); %dy
            dfdt(4) = ki_v_dc*y; %dInt_id

            dfdt(5) = diqdt;  %d Iq_f



    end

    dfdt = dfdt.';
    end