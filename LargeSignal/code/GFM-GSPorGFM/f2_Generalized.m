function dfdt = f2_Generalized(x)

    a1 = evalin('base','a1');
    b1 = evalin('base','b1');
    c1 = evalin('base','c1');
    d1 = evalin('base','d1');
    a2 = evalin('base','a2');
    b2 = evalin('base','b2');
    c2 = evalin('base','c2');
    d2 = evalin('base','d2');

   
    delta1 = x(1);
    delta2 = x(2);

 
    dfdt(1) = a1 + b1 * cos(delta2-delta1) + c1 * sin(delta1) + d1 * sin(delta2-delta1);
    dfdt(2) = a2 + b2 * cos(delta1-delta2) + c2 * sin(delta2) + d2 * sin(delta1-delta2);

    dfdt = dfdt.';

    end