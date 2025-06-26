function dfdt = f_DVCslow(x)


    Iq = evalin('base','Iq');  %for constant Iq
    Ug =evalin('base','Ug');
        
    kp_v_dc = evalin('base','kp_v_dc');
    ki_v_dc = evalin('base','ki_v_dc');

    Ws = evalin('base','Ws');

    Xg = evalin('base','Xg');
    Rg = evalin('base','Rg');
    C_dc = evalin('base','C_dc');
    Pref = evalin('base','Pref');

    Id = x(1);
    y = x(2); %vdc^2-vdc_ref^2


    sindelta = (Id*Xg+Iq*Rg)/Ug;
    cosdelta = sqrt(1-sindelta^2);
    if (Id<Ug/Xg) && (Id>-Ug/Xg)
        P = Id*Ug*cosdelta - Iq*Ug*sindelta + (Id^2+Iq^2)*Rg;
    else
        P=0;
    end

    dfdt(1) = ki_v_dc*y+kp_v_dc*2/(C_dc)*(Pref-P);
    dfdt(2) = 2/(C_dc)*(Pref-P);

    dfdt = dfdt.';

    end