function dfdt = f2_GFMGFL(x)

    model = evalin('base','model');
    I1d = evalin('base','I1d');
    I1q = evalin('base','I1q');
    Ug =evalin('base','Ug');

    k1p = evalin('base','k1p');

    Ws = 50*2*pi;

    X1 = evalin('base','X1');
    X2 = evalin('base','X2');
    Xg = evalin('base','Xg');


    Vgfm = evalin('base','Vgfm');
    Kgfm = evalin('base','Kgfm');
    Pref = evalin('base','Pref');

    X2pg = X2*Xg/(X2+Xg);
    L2pg = X2pg/Ws;
    L1 = X1/Ws;



   
    delta1 = x(1);
    delta2 = x(2);
    

    R_ratio = 0.1;




    if model == "precise"
        dfdt(1) = (k1p/(1-k1p*(L1+L2pg)*I1d)) * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1));
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(Xg+X2) + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));


    elseif model =="rough"   %"rough" ignore w2 & w1
        dfdt(1) = k1p * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1));
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(Xg+X2) + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));

    else %resistance 
        dfdt(1) = k1p * ((X1+X2pg)*I1d-Xg/(X2+Xg)*Vgfm*sin(delta1-delta2)-X2/(X2+Xg)*Ug*sin(delta1)-I1q*R_ratio*(X1+X2pg));
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)*(Xg+X2)/((Xg+X2)*(1+R_ratio))^2 -(Vgfm^2-Vgfm*Ug*cos(delta2))*(Xg+X2)*R_ratio/((Xg+X2)*(1+R_ratio))^2 + Xg/(Xg+X2)*I1d*Vgfm*cos(delta1-delta2) - Xg/(Xg+X2)*I1q*Vgfm*sin(delta1-delta2));
    end
    dfdt = dfdt.';


end





