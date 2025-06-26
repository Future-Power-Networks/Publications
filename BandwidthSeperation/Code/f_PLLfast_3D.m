function dfdt = f_PLLfast_3D(x)

    Ug =evalin('base','Ug');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');

    Id = x(1);
    y = x(2); %vdc^2-vdc_ref^2
    Iq = x(3);

    

    sindelta = (Id*Xg+Iq*Rg)/Ug;
    if sindelta >=1
        sindelta = 1;
        cosdelta = 0;
    elseif sindelta <= -1
        sindelta = -1;
        cosdelta = 0;
    else
        cosdelta = sqrt(1-sindelta^2);
    end


    
    Vd = Ug*cosdelta - Iq*Xg;


    P = Id*Ug*cosdelta - Iq*Ug*sindelta;
    
    w_tvc = evalin('base','w_tvc');
    m_v = evalin('base','m_v');
    Vac_ref = evalin('base','Vac_ref');



    dfdt(1) = ki_v_dc*y+kp_v_dc*2/(C_dc)*(Pref-P);
    dfdt(2) = 2/(C_dc)*(Pref-P);
    dfdt(3) = w_tvc*(-Iq + m_v*(Vd-Vac_ref));

    dfdt = dfdt.';

end