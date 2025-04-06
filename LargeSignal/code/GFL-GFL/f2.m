function dfdt = f2(x)

    V = 1;

    X10 = 1;
    X20 = 1;
    X12 = 0.5;

    P1 = 0.8;
    P2 = 0.8;

%     X10 = 1;
%     X20 = 1;
%     X12 = 10;
% 
%     P1 = 0.8;
%     P2 = 0.8;

    D1 = 100*2*pi;
    D2 = 1*2*pi;

    K1 = V^2 / X10;
    K2 = V^2 / X20;
    K12 = V^2 / X12;
    K21 = K12;
    
    delta1 = x(1);
    delta2 = x(2);
    
    dfdt(1) = (D1^-1) * (P1 - K1*sin(delta1) - K12*sin(delta1-delta2));
    dfdt(2) = (D2^-1) * (P2 - K2*sin(delta2) - K21*sin(delta2-delta1));

    dfdt = dfdt.';
end