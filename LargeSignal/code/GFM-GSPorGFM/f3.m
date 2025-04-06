function dfdt = f3(x)

    V = 1;
    X10 = 1;
    X20 = 1;
    X30 = 1;

    X12 = 5;
    X23 = 5;
    X31 = 5;
    
    K1 = V^2 / X10;
    K2 = V^2 / X20;
    K3 = V^2 / X30;
    
    K12 = V^2 / X12;
    K23 = V^2 / X23;
    K31 = V^2 / X31;
    K21 = K12;
    K32 = K23;
    K13 = K31;

    P1 = 0.8;
    P2 = 0.8;
    P3 = 0.0;
    
    D1 = 100;
    D2 = 100;
    D3 = 100;
    
    delta1 = x(1);
    delta2 = x(2);
    delta3 = x(3);
    
    dfdt(1) = (D1^-1) * (P1 - K1*sin(delta1) - K12*sin(delta1-delta2) - K13*sin(delta1-delta3));
    dfdt(2) = (D2^-1) * (P2 - K2*sin(delta2) - K21*sin(delta2-delta1) - K23*sin(delta2-delta3));
    dfdt(3) = (D3^-1) * (P3 - K3*sin(delta3) - K31*sin(delta3-delta1) - K32*sin(delta3-delta2));

    dfdt = dfdt.';
end