function [V,Vd]=voltage_cal(delta1,delta2,I2q)
    I1d = evalin('base','I1d');
    I1q = evalin('base','I1q');
    Ug =evalin('base','Ug');
    Xg = evalin('base','Xg');
    X2 = evalin('base','X2');
    X2sum = Xg + X2;
    
    Vd=Ug*cos(delta2)-X2sum*I2q-Xg*I1d*sin(delta1-delta2)-Xg*I1q*cos(delta1-delta2);
    Vq=-Ug*sin(delta2)+Xg*I1d*cos(delta1-delta2)-Xg*I1q*sin(delta1-delta2);
    V=sqrt(Vd^2+Vq^2);
end