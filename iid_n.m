clear; clc

%define parameters

n = 8; % # nodes
T = 50; % # periods forward running the risk free rate
%m=12; %number of edges, maximum n(n-1)
risk_ratio = 1;
m = .015;
lambda= m; %internsity of distress
eta=m*risk_ratio; %probability of recovery
delta = .95; %time prefrence
gamma = 1; %RRA parameter of CRRA

%Dividend at each state
x(1) = 1; %x_lowerbar, value of idiosyncratic component at busat state
x(2) = 1e10; %x_bar, value of idiosyncratic component at boom state

Delta = zeros(n); %complete network, non-directional
%%
H = zeros(n,2^n); %States of the world (H)
Lambda = zeros(n,2^n); %collect all lambda_i for each node at each state
Pi = zeros(2^n,2^n); %markov matrix
D = NaN(2^n,1); %total dividend at each state

for k=1:2^n %state of the world
    str = dec2bin(k-1,n); %generate binary number k-1
    digit = str-'0'; %arrange each digit in a row vector
    H(:,k) = digit'; %define this vector as state of the world k
    D(k) = x(1)*sum(H(:,k))+x(2)*(n-sum(H(:,k))); %total dividend at k

    %Lambda(:,k) = lambda*Delta*H(:,k)+lambda; %With natural distress rate
    Lambda(:,k) = lambda*Delta*H(:,k); %Without natural rate
end
%%

for k_origin = 1:2^n
    for k_dest = 1:2^n

        aux = NaN(n,1);
        for i = 1:n
aux(i) = (1-eta)*H(i,k_origin)*H(i,k_dest)+eta*H(i,k_origin)*(1-H(i,k_dest))+...
         Lambda(i,k_origin)*(1-H(i,k_origin))*H(i,k_dest)+(1-Lambda(i,k_origin))*(1-H(i,k_origin))*(1-H(i,k_dest));
        end
        Pi(k_origin,k_dest) = prod(aux);
        
    end
end


MU = @(u) u.^(-gamma);
A = Pi.*(MU(D)*(1./MU(D))')'*delta;
%MU = @(u) exp(-delta)*u.^(-gamma);
%A = Pi.*(MU(D)*(1./MU(D))');

r = NaN(T,2^n);
for t = 1:T
    r(t,:) = (1./(sum(mpower(A,t),2))).^(1/t);
    %r(t,:) = (1./(sum(mpower(A,t),2))-1)/t+1;
    %r(t,:) = (-log(sum(mpower(A,t),2)))/t+1;
end

%%

%{
**************************************************************************
                                FIGURES
**************************************************************************
%}

iid_fig = figure;
leg = cell(1,n+1);
%Center not distressed, i-1 periphery distressed

% there are four colors, so grouping n states to 4 groups we have x in each
% group: x = ceil((n-2)/4)
state_num_in_group = ceil((n-2)/4);
for i = 1:n+1
    if i == 1
        color_used = 'g--';
    elseif i > 1 && i <= state_num_in_group+1
        color_used = 'c-';
    elseif i > 1*state_num_in_group+1 && i <= 2*state_num_in_group+1
        color_used = 'b-';
    elseif i > 2*state_num_in_group+1 && i <= 3*state_num_in_group+1
        color_used = 'k-';
    elseif i > 3*state_num_in_group+1 && i <= 4*state_num_in_group+1 && ...
            i ~= n+1
        color_used = 'm-';
    elseif i == n+1
        color_used = 'r-';        
    end
    
    leg{i} = num2str(i-1);
    if i == 1
        plot(r(:,1)-1,color_used) % 0 distressed        
    else
        plot(r(:,2^(n-1)+2^(i-2))-1,color_used) % i-1 distressed
    end
    hold on
end


title(sprintf('i.i.d. network, n=%d',n))
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))
axis([0,50,0,0.13]) %for comparison!



legend(leg,'Location','northeast','Orientation','vertical')

cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(iid_fig,'-dpng','-r100',...
    sprintf('iid_net_n%d_lambda=%.2f.jpg',...
    n,lambda))













