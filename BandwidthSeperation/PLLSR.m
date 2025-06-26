for Id = 2.05

%%
x=(0:0.1:1)*2*pi;
n = length(x);
x_set = zeros(2,n);
x_set(1,:) = x;

torralence = 1e-2; 
mm = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',100000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            vv = V(:,~(~sig)); 
            
            ep_set(mm).xep = xep; %#ok<*SAGROW> 
            ep_set(mm).A = A;
            ep_set(mm).Lambda = Lambda;
            ep_set(mm).V = V;   
            ep_set(mm).v = v;     % stable eigenvectors of unstable ep 
            ep_set(mm).vv = vv;
            ep_set(mm).flag = flag;
           
            mm = mm+1;
        end
    end
end

for mm = 1:length(ep_set)
    disp_v('Index',mm);
    disp_v('Equilibrium',ep_set(mm).xep);
    disp_v('Eigenvalue', ep_set(mm).Lambda);
    disp_v('Eigenvector',ep_set(mm).V);
end

clear ep_set_ext;
for n = 1:length(ep_set)
    mm = (n-1)*10;
    ep_set_ext(mm+1)=ep_set(n); %#ok<*AGROW> 
    ep_set_ext(mm+2)=ep_set(n);
    ep_set_ext(mm+3)=ep_set(n); 
    ep_set_ext(mm+4)=ep_set(n); 
    ep_set_ext(mm+5)=ep_set(n);
    ep_set_ext(mm+6)=ep_set(n);
    ep_set_ext(mm+7)=ep_set(n);
    ep_set_ext(mm+8)=ep_set(n);
    ep_set_ext(mm+9)=ep_set(n);
    ep_set_ext(mm+10)=ep_set(n);
    ep_set_ext(mm+2).xep(1) = ep_set(n).xep(1) - 2*pi;
    ep_set_ext(mm+3).xep(1) = ep_set(n).xep(1) + 2*pi;
    ep_set_ext(mm+4).xep(1) = ep_set(n).xep(1) + 4*pi;
    ep_set_ext(mm+5).xep(1) = ep_set(n).xep(1) + 6*pi;
    ep_set_ext(mm+6).xep(1) = ep_set(n).xep(1) + 10*pi;
    ep_set_ext(mm+7).xep(1) = ep_set(n).xep(1) + 18*pi;
    ep_set_ext(mm+8).xep(1) = ep_set(n).xep(1) + 26*pi;
    ep_set_ext(mm+9).xep(1) = ep_set(n).xep(1) + 34*pi;
    ep_set_ext(mm+10).xep(1) = ep_set(n).xep(1) - 4*pi;
end
%%
f10=figure(10);
hold on;
grid on;
ymin=-5;
ymax=20;
color_code = {'black','magenta','red','black'};
axis([-3/2*pi,3/2*pi,ymin,ymax]);
xticks(-2*pi:pi/2:2*pi);
xticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 14);
for mm = 1 : length(ep_set_ext)
    xep = ep_set_ext(mm).xep;
    flag= ep_set_ext(mm).flag;
    scatter(xep(1),xep(2),color_code{flag+1},'LineWidth', 1.5);
    if flag == 1
        v = ep_set_ext(mm).v;
        vv = ep_set_ext(mm).vv;
        perturb = 1e-3;
        [~ , x_p] = ode78(@f_backward,[0,80],xep+v*perturb,odeset('RelTol',1e-5));
        [~ , x_n] = ode78(@f_backward,[0,80],xep-v*perturb,odeset('RelTol',1e-5));
        x_all = [flip(x_n,1);x_p];
        plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);
    end
end
%%
% figure(f10);
% plot(delta_fault,Int_fault,'k-','linewidth',1.5);
% plot(delta_post(1),Int_post(1),'k.','MarkerSize',15);
% plot(delta_post,Int_post,'k-','linewidth',1.5);

end




%% function
function yes = isnewxep(ep_set,xep,torr)
    if isempty(ep_set)
        yes = 1;
        return;
    end
    minerr = inf;
    for m = 1 : length(ep_set)
        err = abs(xep - ep_set(m).xep);
        err = min(err, abs(2*pi-err));
        err = max(err);
        if minerr > err
            minerr = err;
        end
    end
    if(minerr>torr)
        yes = 1;
    else
        yes = 0;
    end
end



function dfdt = f(x)
        dfdt = f_GFL(x);
end

function out = maxabs(in)

    out = abs(in);
    
    while length(out) > 1
        out = max(out);
    end

end

function disp_v(msg,v)
    disp([msg '=']);
    disp(v);
end

function dfdt = f_backward(t,x)
    dfdt = -f(x);
end

