function dfdt = f2_GFMGFLV2(x)

    model = evalin('base','model');
    I1d = evalin('base','I1d');
    Ug =evalin('base','Ug');
    k1p = evalin('base','k1p');
    Kq = evalin('base','Kq');

    Ws = 50*2*pi;

    X1 = evalin('base','X1');
    X2 = evalin('base','X2');
    Xg = evalin('base','Xg');
    R1 = evalin('base','R1');
    R2 = evalin('base','R2');
    Rg = evalin('base','Rg');


    Vgfm = evalin('base','Vgfm');
    Kgfm = evalin('base','Kgfm');
    Pref = evalin('base','Pref');

    V1ref = evalin('base','V1ref');

    X2pg = X2*Xg/(X2+Xg);
    if R2+Rg>0
       R2pg = R2*Rg/(R2+Rg);
    else
       R2pg = 0;
    end
    L2pg = X2pg/Ws;
    L1 = X1/Ws;
    XX = Xg*X1+X2*X1+X2*Xg;

    Z1 = R1+X1*1j;
    Z2 = R2+X2*1j;
    Zg = Rg+Xg*1j;
    Z2sum = Z2+Zg;



   
    delta2 = x(1);
    delta1 = x(2);

    %I1q = kq*(Xg/(X2+Xg)*Vgfm*cos(delta1-delta2)+X2/(X2+Xg)*Ug*cos(delta1)-V1ref)/(1+kq*(XX/(X2+Xg)));
    I1q = Kq*( real(Zg/(Z2+Zg)*Vgfm*exp((delta2-delta1)*1j)) + real( Z2/(Z2+Zg)*Ug*exp(-1*1j*delta1) )-V1ref +(R1+R2pg)*I1d)/(1+Kq*(XX/(X2+Xg)));

    if model == "precise"
        dfdt(1) = (k1p/(1-k1p*(L1+L2pg)*I1d)) * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1));
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(Xg+X2) + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));


    elseif model == "rough" %ignore w2 & w1
%         dfdt(1) = k1p * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1));
%         dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(Xg+X2) + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));
        dfdt(2) = k1p * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1));
        dfdt(1) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(Xg+X2) + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));

    else
          dfdt(2) = k1p * ((X1+X2pg)*I1d+imag(Zg/Z2sum*Vgfm*exp((delta2-delta1)*1j))+imag(Z2/Z2sum*Ug*exp(delta1*-1j))+I1q*(R1+R2pg)); 
          dfdt(1) = Kgfm*(Pref +real(Zg/Z2sum*I1d*Vgfm*exp((delta1-delta2)*1j)) + real(Zg/Z2sum*I1q*Vgfm*exp((delta1-delta2+pi/2)*1j)) - real(Vgfm^2/Z2sum-Vgfm*Ug/Z2sum*exp((-delta2)*1j)) );

    end
    dfdt = dfdt.';


end





