function dfdt = f2_GFL(x)

    model = evalin('base','model');
    I1d = evalin('base','I1d');
    I1q = evalin('base','I1q');
    I2d = evalin('base','I2d');
    I2q = evalin('base','I2q');
    Ug =evalin('base','Ug');

    k1p = evalin('base','k1p');
    k2p = evalin('base','k2p');

    Ws = 50*2*pi;

    X1 = evalin('base','X1');%X10*X12/(X10+X20+X12);
    X2 = evalin('base','X2');%X20*X12/(X10+X20+X12);
    Xg = evalin('base','Xg');%X10*X20/(X10+X20+X12);

    w_limit =evalin('base','w_limit');



    X1sum = X1+Xg;
    X2sum = X2+Xg;
    L1sum = X1sum/Ws;
    L2sum = X2sum/Ws;
    Lg = Xg/Ws;


    M1 = (1-k1p*L1sum*I1d)/k1p;
    M2 = (1-k2p*L2sum*I2d)/k2p;
   
    delta1 = x(1);
    delta2 = x(2);

    C1 = I2d*cos(delta2-delta1)-I2q*sin(delta2-delta1);
    C2 = I1d*cos(delta1-delta2)-I1q*sin(delta1-delta2);
    MM1 = M1 -1/M2*Lg^2*C1*C2;
    MM2 = M2 -1/M1*Lg^2*C2*C1;

    if model == "precise"
    
            dfdt(1) = (MM1^-1) * (X1sum*I1d + Xg*C1-Ug*sin(delta1) + Lg*C1/M2*(X2sum*I2d + Xg*C2-Ug*sin(delta2)) );
            dfdt(2) = (MM2^-1) * (X2sum*I2d + Xg*C2-Ug*sin(delta2) + Lg*C2/M1*(X1sum*I1d + Xg*C1-Ug*sin(delta1)) );

    elseif model == "approximation1" %regards w1 = w2
            
            dfdt(1) = ((M1-Lg*C1)^-1) * (X1sum*I1d + Xg*C1-Ug*sin(delta1));
            dfdt(2) = ((M2-Lg*C2)^-1) * (X2sum*I2d + Xg*C2-Ug*sin(delta2));

    elseif model == "approximation2" % ignore w2

            dfdt(1) = (M1^-1) * (X1sum*I1d + Xg*C1-Ug*sin(delta1));
            dfdt(2) = (M2^-1) * (X2sum*I2d + Xg*C2-Ug*sin(delta2));

    elseif model == "rough" % ignore w2 & w1

            dfdt(1) = k1p * ( X1sum*I1d + Xg*C1-Ug*sin(delta1));%
            dfdt(2) = k2p * ( X2sum*I2d + Xg*C2-Ug*sin(delta2));%
            % omega-limit
        if dfdt(1)>w_limit  %not correct
            dfdt(1)=w_limit;
        elseif dfdt(1)<-w_limit
            dfdt(1)=-w_limit;
        end
        if dfdt(2)>w_limit
            dfdt(2)=w_limit;
        elseif dfdt(2)<-w_limit
            dfdt(2)=-w_limit;
        end

    elseif model == "precise_limit" % ignore w2 & w1   
        dfdt(1) = (MM1^-1) * (X1sum*I1d + Xg*C1-Ug*sin(delta1) + Lg*C1/M2*(X2sum*I2d + Xg*C2-Ug*sin(delta2)) );
        dfdt(2) = (MM2^-1) * (X2sum*I2d + Xg*C2-Ug*sin(delta2) + Lg*C2/M1*(X1sum*I1d + Xg*C1-Ug*sin(delta1)) );
        if dfdt(1)>w_limit
            dfdt(1)=w_limit;
        elseif dfdt(1)<-w_limit
            dfdt(1)=-w_limit;
        end
% 
        if dfdt(2)>w_limit
            dfdt(2)=w_limit;
        elseif dfdt(2)<-w_limit
            dfdt(2)=-w_limit;
        end
    
    else % frequency limit & rough
        dfdt(1) = k1p * (X1sum*I1d + Xg*C1-Ug*sin(delta1));
        dfdt(2) = k2p * (X2sum*I2d + Xg*C2-Ug*sin(delta2));
% 
        if dfdt(1)>w_limit  %not correct
            dfdt(1)=w_limit;
        elseif dfdt(1)<-w_limit
            dfdt(1)=-w_limit;
        end
% 
        if dfdt(2)>w_limit
            dfdt(2)=w_limit;
        elseif dfdt(2)<-w_limit
            dfdt(2)=-w_limit;
        end
    end
    dfdt = dfdt.';

    end