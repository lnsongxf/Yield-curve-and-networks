clear; clc



%define parameters

n = 9; % # nodes
T = 50; % # periods forward running the risk free rate
m=12; %number of edges, maximum n(n-1)
lambda= .06; %internsity of distress
eta=.02; %probability of recovery

delta = .95; %time prefrence
gamma = 2; %RRA parameter of CRRA

%Dividend at each state
x(1) = 1; %x_lowerbar, value of idiosyncratic component at busat state
x(2) = 2; %x_bar, value of idiosyncratic component at boom state

Delta = ones(n)-eye(n); %complete network, non-directional

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
A = Pi.*(MU(D)*(1./MU(D))')*delta;
%MU = @(u) exp(-delta)*u.^(-gamma);
%A = Pi.*(MU(D)*(1./MU(D))');

r = NaN(T,2^n);
for t = 1:T
    r(t,:) = (1./(sum(mpower(A,t),2))).^(1/t);
end

plot(r-1)
%%
%Figure adapted for n=4 case
cmp_fig = figure;
plot(r(:,1)-1,'k') % 0 distressed
hold on
plot(r(:,9)-1,'b--') % 1 distressed
hold on
plot(r(:,10)-1,'r-.') % 2 distressed
hold on
plot(r(:,12)-1,'g.-') % 3 distressed
hold on
plot(r(:,16)-1,'c.-.') % 4 distressed

title(sprintf('Complete network, n=%d',n))
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))
axis([0,50,0,0.13])

legend('0 distressed','1 distressed','2 distressed',...
    '3 distressed','4 distressed',...
    'Location','northeast','Orientation','vertical')
%}    
%%
%{
cd C:\Users\Oren\Documents\MATLAB\Network
print(cmp_fig,'-djpeg',...
    sprintf('cmp-eta%.2f-lambda%.2f',eta,lambda))
%}







