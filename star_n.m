clear; clc



%define parameters

n = 8; % # nodes
T = 50; % # periods forward running the risk free rate
m=12; %number of edges, maximum n(n-1)
lambda= .06; %internsity of distress
eta=.02; %probability of recovery

delta = .95; %time prefrence
gamma = 2; %RRA parameter of CRRA

%Dividend at each state
x(1) = 1; %x_lowerbar, value of idiosyncratic component at busat state
x(2) = 2; %x_bar, value of idiosyncratic component at boom state

Delta = star_net(n,1); %star network, two way eff
%Delta = [0,1,1,1;zeros(3,4)]'; %star network, one way eff
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

    %Lambda(:,k) = lambda*Delta*H(:,k)+lambda; %With natural prob
    Lambda(:,k) = lambda*Delta*H(:,k); %Without natural prob
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
end


%%

%Colors = {'g','b','c','k','r'};
star_fig = figure;
for i = 1:n+1
    if i == 1
        color_used = 'g';
    elseif i == 2 || i == 3
        color_used = 'b';
    elseif i == 4 || i == 5
        color_used = 'c';
    elseif i == 6 || i == 7
        color_used = 'k';
    else
        color_used = 'r';
    end
    plot(r(:,2^(i-1))-1,color_used) % i-1 periphery firm is distressed; center not
    hold on
end

for i = 1:n
    
    if i == 1 || i == 2
        color_used = 'b--';
    elseif i == 3 || i == 4
        color_used = 'c--';
    elseif i == 5 || i == 6
        color_used = 'k--';
    elseif i == 7 || i == 8
        color_used = 'r--';
    end
    
    plot(r(:,2^(n-1)+2^(i-1))-1,color_used) % 0 periphery firm is distressed; center yes
    hold on
end
%{
plot(r(:,2^0)-1,'g') % 0 periphery firm is distressed; center not
hold on
plot(r(:,2^1)-1,'b') % 1 periphery firm is distressed; center not
hold on
plot(r(:,2^2)-1,'b') % 2 periphery firm is distressed; center not
hold on
plot(r(:,2^3)-1,'c') % 3 periphery firm is distressed; center not
hold on
plot(r(:,2^4)-1,'c') % 4 periphery firm is distressed; center not
hold on
plot(r(:,2^5)-1,'k') % 5 periphery firm is distressed; center not
hold on
plot(r(:,2^6)-1,'k') % 6 periphery firm is distressed; center not
hold on
plot(r(:,2^7)-1,'r') % 7 periphery firm is distressed; center not
hold on

plot(r(:,2^(n-1)+2^1)-1,'b--') % 1 periphery firm is distressed; center yes
hold on
plot(r(:,2^(n-1)+2^2)-1,'b--') % 2 periphery firm is distressed; center yes
hold on
plot(r(:,2^(n-1)+2^3)-1,'c--') % 3 periphery firm is distressed; center yes
hold on
plot(r(:,2^(n-1)+2^4)-1,'c--') % 4 periphery firm is distressed; center yes
hold on
plot(r(:,2^(n-1)+2^5)-1,'k--') % 5 periphery firm is distressed; center yes
hold on
plot(r(:,2^(n-1)+2^6)-1,'k--') % 6 periphery firm is distressed; center yes
hold on
plot(r(:,2^(n-1)+2^7)-1,'r--') % 7 periphery firm is distressed; center yes
hold on
%}
title(sprintf('Star network, n=%d',n))
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))
%axis([0,50,0,0.13])


leg = cell(1,2*n);

for i = 1:n
    leg{i} = strcat('0',num2str(i-1));
    leg{i+n} = strcat('1',num2str(i-1));
end

legend(leg,'Location','northeast','Orientation','vertical')

cd C:\Users\Oren\Documents\MATLAB\Network\figures
print(star_fig,'-dpng','-r100',...
    sprintf('star_net_n%d_lambda=%.2f.jpg',...
    n,lambda))


%{
legend('No distress','1 periphery distressed',...
    '2 periphery distressed','3 periphery distressed',...
    'Only center distressed','Center & 1 periphery',...
    'Center & 2 periphery','Center & 3 periphery',...
    'Location','northeast','Orientation','vertical')
%}    
%%









