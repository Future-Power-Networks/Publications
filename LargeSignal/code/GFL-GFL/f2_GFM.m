function dfdt = f2_GFM(x)
    model = evalin('base','model');

    Vgfm1 = evalin('base','Vgfm1');
    Vgfm2 = evalin('base','Vgfm2');
    Kgfm1 = evalin('base','Kgfm1');
    Kgfm2 = evalin('base','Kgfm2');
    Pref1 =evalin('base','Pref1');
    Pref2 =evalin('base','Pref2');
    Ug = evalin('base','Ug');

    Ws = 50*2*pi;


    

    X12 = evalin('base','X12');
    X1g = evalin('base','X1g');
    X2g = evalin('base','X2g');


    R_ratio = 0.1;
    
    delta1 = x(1);
    delta2 = x(2);
    if model == "resistance"
        dfdt(1) = Kgfm1 * (Pref1 - Vgfm1*Vgfm2*X12/(X12*(1+R_ratio))^2*sin(delta1-delta2) - (Vgfm1^2-Vgfm1*Vgfm2*cos(delta1-delta2))*X12*R_ratio/(X12*(1+R_ratio))^2 - Vgfm1*Ug*sin(delta1)*X1g/(X1g*(1+R_ratio))^2 -(Vgfm1^2-Vgfm1*Ug*cos(delta1))*X1g*R_ratio/(X1g*(1+R_ratio))^2);
        dfdt(2) = Kgfm2 * (Pref2 - Vgfm1*Vgfm2*X12/(X12*(1+R_ratio))^2*sin(delta2-delta1) - (Vgfm2^2-Vgfm1*Vgfm2*cos(delta2-delta1))*X12*R_ratio/(X12*(1+R_ratio))^2 - Vgfm2*Ug*sin(delta2)*X2g/(X2g*(1+R_ratio))^2 -(Vgfm2^2-Vgfm2*Ug*cos(delta2))*X2g*R_ratio/(X2g*(1+R_ratio))^2);
    else
        dfdt(1) = Kgfm1 * (Pref1 - Vgfm1*Vgfm2/X12*sin(delta1-delta2) - Vgfm1*Ug/X1g*sin(delta1));
        dfdt(2) = Kgfm2 * (Pref2 - Vgfm1*Vgfm2/X12*sin(delta2-delta1) - Vgfm2*Ug/X2g*sin(delta2));
    end

    dfdt = dfdt.';
end
