function dfdt = f2_GFM_fault(x)
    model = evalin('base','model');

    Vgfm1 = evalin('base','Vgfm1');
    Vgfm2 = evalin('base','Vgfm2');
    Kgfm1 = evalin('base','Kgfm1');
    Kgfm2 = evalin('base','Kgfm2');
    Pref1 =evalin('base','Pref1');
    Pref2 =evalin('base','Pref2');
    Ug = evalin('base','Ug_fault');

    Ws = 50*2*pi;



    X1 = evalin('base','X1f');
    X2 = evalin('base','X2f');
    Xg = evalin('base','Xg_f');
    R1 = evalin('base','R1f');
    R2 = evalin('base','R2f');
    Rg = evalin('base','Rg_f');

    Z1 = R1+X1*1j;
    Z2 = R2+X2*1j;
    Zg = Rg+Xg*1j;
    ZZ = Z1*Z2+Z2*Zg+Z1*Zg;
    Z12 = ZZ/Zg;
    Z1g = ZZ/Z2;
    Z2g = ZZ/Z1;
    Z11 = ZZ/(Z1+Zg);
    Z22 = ZZ/(Z1+Zg);
    X12 = (X1*X2+X2*Xg+X1*Xg)/Xg;
    X1g = (X1*X2+X2*Xg+X1*Xg)/X2;
    X2g = (X1*X2+X2*Xg+X1*Xg)/X1;

    
    delta1 = x(1);
    delta2 = x(2);
    if model == "rough"
%         dfdt(1) = Kgfm1 * (Pref1 - Vgfm1*Vgfm2*X12/(X12*(1+R_ratio))^2*sin(delta1-delta2) - (Vgfm1^2-Vgfm1*Vgfm2*cos(delta1-delta2))*X12*R_ratio/(X12*(1+R_ratio))^2 - Vgfm1*Ug*sin(delta1)*X1g/(X1g*(1+R_ratio))^2 -(Vgfm1^2-Vgfm1*Ug*cos(delta1))*X1g*R_ratio/(X1g*(1+R_ratio))^2);
%         dfdt(2) = Kgfm2 * (Pref2 - Vgfm1*Vgfm2*X12/(X12*(1+R_ratio))^2*sin(delta2-delta1) - (Vgfm2^2-Vgfm1*Vgfm2*cos(delta2-delta1))*X12*R_ratio/(X12*(1+R_ratio))^2 - Vgfm2*Ug*sin(delta2)*X2g/(X2g*(1+R_ratio))^2 -(Vgfm2^2-Vgfm2*Ug*cos(delta2))*X2g*R_ratio/(X2g*(1+R_ratio))^2);
        dfdt(1) = Kgfm1 * (Pref1 - Vgfm1*Vgfm2/X12*sin(delta1-delta2) - Vgfm1*Ug/X1g*sin(delta1));
        dfdt(2) = Kgfm2 * (Pref2 - Vgfm1*Vgfm2/X12*sin(delta2-delta1) - Vgfm2*Ug/X2g*sin(delta2));
    else %resistance
        dfdt(1) = Kgfm1 * (Pref1 - real(Vgfm1^2/Z11-Vgfm1*Ug/Z1g*exp(-1*1j*delta1)- Vgfm1*Vgfm2/Z12*exp(1j*(delta2-delta1))));
        dfdt(2) = Kgfm2 * (Pref2 - real(Vgfm2^2/Z22-Vgfm2*Ug/Z2g*exp(-1*1j*delta2)- Vgfm1*Vgfm2/Z12*exp(1j*(delta1-delta2))));
    end

    dfdt = dfdt.';
end
