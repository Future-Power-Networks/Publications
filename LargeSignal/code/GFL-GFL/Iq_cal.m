function I2q=Iq_cal(delta1,delta2)

    I1d = evalin('base','I1d');
    I1q = evalin('base','I1q');
    Ug =evalin('base','Ug');
    Xg = evalin('base','Xg');
    X2 = evalin('base','X2');
    X2sum = Xg + X2;

    V1ref = evalin('base','V1ref');
    kq = evalin('base','Kq');
    I2q = kq*(Ug*cos(delta2)-V1ref-Xg*I1q*cos(delta1-delta2)-Xg*I1d*sin(delta1-delta2))/(1+kq*X2sum);

end