function V = voltage1(x)
    I1d = evalin('base','I1d');
    I1q = evalin('base','I1q');
    I2d = evalin('base','I2d');
    I2q = evalin('base','I2q');
    Ug =evalin('base','Ug');

    Ws = 50*2*pi;

    X1 = evalin('base','X1');
    X2 = evalin('base','X2');
    Xg = evalin('base','Xg');

    ph1 = atan(I1q/I1d);
    ph2 = atan(I2q/I2d);
    delta1 = x(1);
    delta2 = x(2);
    I1 = sqrt(I1d^2+I1q^2);
    I2 = sqrt(I2d^2+I2q^2);


    V_square = (X1+Xg)^2*I1^2+Ug^2+Xg^2*I2^2-2*Ug*(X1+Xg)*I1*sin(ph1+delta1)-2*Ug*Xg*I2*sin(ph2+delta2)+2*(X1+Xg)*I1*Xg*I2*cos(ph1+delta1-ph2-delta2);
    V = sqrt(V_square);

end