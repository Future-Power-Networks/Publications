function dfdt = f2_STATCOM(x)

    model = evalin('base','model');
    I1d = evalin('base','I1d');
    I1q = evalin('base','I1q');
    I2d = evalin('base','I2d');
    Ug =evalin('base','Ug');

    k1p = evalin('base','k1p');
    k2p = evalin('base','k2p');

    Ws = 50*2*pi;

    X1 = evalin('base','X1');%X10*X12/(X10+X20+X12);
    X2 = evalin('base','X2');%X20*X12/(X10+X20+X12);
    Xg = evalin('base','Xg');%X10*X20/(X10+X20+X12);

    w_limit =evalin('base','w_limit');

    V1ref = evalin('base','V1ref');
    kq = evalin('base','Kq');
    




    X1sum = X1+Xg;
    X2sum = X2+Xg;
    L1sum = X1sum/Ws;
    L2sum = X2sum/Ws;
    Lg = Xg/Ws;


    M1 = (1-k1p*L1sum*I1d)/k1p;
    M2 = (1-k2p*L2sum*I2d)/k2p;
   
    delta1 = x(1);
    delta2 = x(2);




    if model == "precise"
            I2q = kq*(Ug*cos(delta2)-V1ref-Xg*I1q*cos(delta1-delta2)-Xg*I1d*sin(delta1-delta2))/(1+kq*X2sum);
            C1 = I2d*cos(delta2-delta1)-I2q*sin(delta2-delta1);
            C2 = I1d*cos(delta1-delta2)-I1q*sin(delta1-delta2);
            MM1 = M1 -1/M2*Lg^2*C1*C2;
            MM2 = M2 -1/M1*Lg^2*C2*C1;
    
            dfdt(1) = (MM1^-1) * (X1sum*I1d + Xg*C1-Ug*sin(delta1) + Lg*C1/M2*(X2sum*I2d + Xg*C2-Ug*sin(delta2)) );
            dfdt(2) = (MM2^-1) * (X2sum*I2d + Xg*C2-Ug*sin(delta2) + Lg*C2/M1*(X1sum*I1d + Xg*C1-Ug*sin(delta1)) );

    elseif model == "approximation1" %consider w1 & w2
            I2q = kq*(Ug*cos(delta2)-V1ref-Xg*I1q*cos(delta1-delta2)-Xg*I1d*sin(delta1-delta2))/(1+kq*X2sum);
            C1 = I2d*cos(delta2-delta1)-I2q*sin(delta2-delta1);
            C2 = I1d*cos(delta1-delta2)-I1q*sin(delta1-delta2);

            dfdt(1) = k1p * (X1sum*I1d + Xg*C1-Ug*sin(delta1));
            dfdt(2) = k2p * (X2sum*I2d + Xg*C2-Ug*sin(delta2));

    else   %"rough" ignore w2 & w1
            I2q = kq*(Ug*cos(delta2)-V1ref-Xg*I1q*cos(delta1-delta2)-Xg*I1d*sin(delta1-delta2))/(1+kq*X2sum);
            C1 = I2d*cos(delta2-delta1)-I2q*sin(delta2-delta1);
            C2 = I1d*cos(delta1-delta2)-I1q*sin(delta1-delta2);
            dfdt(1) = k1p * (X1sum*I1d + Xg*C1-Ug*sin(delta1));
            dfdt(2) = k2p * (X2sum*I2d + Xg*C2-Ug*sin(delta2));

    end
    dfdt = dfdt.';


end





