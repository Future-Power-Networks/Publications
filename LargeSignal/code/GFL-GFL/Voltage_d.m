function [Vd1, Vd2] = Voltage_d(x)
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

    Vd1 = Ug*cos(delta1)-(X1+Xg)*I1q-Xg*I2*sin(delta2-delta1+ph2);
    Vd2 = Ug*cos(delta2)-(X2+Xg)*I2q-Xg*I1*sin(delta1-delta2+ph1);

end