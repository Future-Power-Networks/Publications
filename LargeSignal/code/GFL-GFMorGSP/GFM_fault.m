function dx = f_GFM(x)

    nodes = length(x)+1;

    switch nodes 
        case 3
            Y = [0   1.2   1/2;...
                 1.2   0   1.2;...
                 1/2   1.2  0];
        
            P  = [-1.6; 0.8 ; 0.8];
        
            D =  [100; 1000; 100]./(2*pi);
        otherwise
            error('Not a valid dimension.');
    end

    V = 1;
    K = V^2 .* Y;

    delta=zeros(nodes,1);
    ddelta=zeros(nodes,1);
    dx=zeros(nodes-1,1);

    delta(1) = 0;
    for n = 1:(nodes-1)
        delta(n+1) = x(n) + delta(1); % delta(n);
    end
    
    for n = 1:nodes
        Pe = 0;
        for m = 1:nodes
            Pe = Pe + K(n,m) * sin(delta(n)-delta(m));
        end
        ddelta(n) = (D(n)^-1) * (P(n) - Pe);
    end

    for n = 1:(nodes-1)
        dx(n) = ddelta(n+1);% - ddelta(n);
    end
end