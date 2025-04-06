function dfdt = f2_GFMGFLV_fault(x)

    model = evalin('base','model');
    I1d = evalin('base','I1d');
    Ug =evalin('base','Ug_fault');

    k1p = evalin('base','k1p');

    kq = evalin('base','Kq');

    Ws = 50*2*pi;

    X1 = evalin('base','X1');
    X2 = evalin('base','X2');
    Xg = evalin('base','Xg');


    Vgfm = evalin('base','Vgfm');
    Kgfm = evalin('base','Kgfm');
    Pref = evalin('base','Pref');

    V1ref = evalin('base','V1ref');

    X2pg = X2*Xg/(X2+Xg);
    L2pg = X2pg/Ws;
    L1 = X1/Ws;
    XX = Xg*X1+X2*X1+X2*Xg;



   
    delta1 = x(1);
    delta2 = x(2);

    I1q = kq*(Xg/(X2+Xg)*Vgfm*cos(delta1-delta2)+X2/(X2+Xg)*Ug*cos(delta1)-V1ref)/(1+kq*(XX/(X2+Xg)));




    if model == "precise"
        dfdt(1) = (k1p/(1-k1p*(L1+L2pg)*I1d)) * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1));
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(Xg+X2) + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));


    else   %"rough" ignore w2 & w1
        dfdt(1) = k1p * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1));
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(Xg+X2) + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));

    end
    dfdt = dfdt.';


end





