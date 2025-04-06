function dfdt = f2_GFMGFLP2(x)

    model = evalin('base','model');
    I1q = evalin('base','I1q');
    Ug =evalin('base','Ug');

    k1p = evalin('base','k1p');

    Ws = 50*2*pi;

    X1 = evalin('base','X1');
    X2 = evalin('base','X2');
    Xg = evalin('base','Xg');


    Vgfm = evalin('base','Vgfm');
    Kgfm = evalin('base','Kgfm');
    Pref = evalin('base','Pref');%GFM
    Pgfl_ref = evalin('base','Pgfl');%GFL

    Kp_p = evalin('base','Kp_p');
    kp_i = evalin('base','Ki_p');

    X2pg = X2*Xg/(X2+Xg);
    L2pg = X2pg/Ws;
    L1 = X1/Ws;
    kg = Xg/(X2+Xg);
    k2 = X2/(X2+Xg);
    XX = Xg*X1+X2*X1+X2*Xg;

   
    i_d = x(1);
    delta2 = x(2);
    A = sqrt((Xg*Vgfm)^2+(X2*Ug)^2+2*Xg*Vgfm*X2*Ug*cos(delta2));
    if i_d >= A/XX
        i_d = A/XX;
    end
    Pgfl = i_d/(Xg+X2)*real(sqrt(A^2-i_d^2*XX^2));
    dPgfldI = 1/(Xg+X2)*(A^2-2*i_d^2*XX^2)/real(sqrt(A^2-i_d^2*XX^2));


    if model == "precise"
        dfdt(1) = kp_i*(Pgfl_ref-Pgfl);%/(1+Kp_p*dPgfldI);
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(X2+Xg)+i_d/(X2+Xg)/A^2*(Xg^2*Vgfm^2*sqrt(A^2-XX^2*i_d^2)+X2*Ug*Xg*Vgfm*sin(delta2)*i_d*XX));


    else   %"rough" ignore w2 & w1
        dfdt(1) = kp_i*(Pgfl_ref-Pgfl);%/(1+Kp_p*dPgfldI);
        dfdt(2) = Kgfm*(Pref-Vgfm*Ug*sin(delta2)/(X2+Xg)+i_d/(X2+Xg)/A^2*(Xg^2*Vgfm^2*sqrt(A^2-XX^2*i_d^2)+X2*Ug*Xg*Vgfm*sin(delta2)*i_d*XX));

    end
    dfdt = dfdt.';


end





