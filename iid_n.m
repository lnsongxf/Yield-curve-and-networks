clear; clc

%define parameters

n = 8; % # nodes
T = 50; % # periods forward running the risk free rate
%m=12; %number of edges, maximum n(n-1)
risk_ratio = .95;
m = .4;
lambda= m; %internsity of distress
eta=m*risk_ratio; %probability of recovery

delta = .95; %time prefrence
gamma = 4; %RRA parameter of CRRA

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

%need to find n+1 states out of 2^n according to the sum of firms being
%distressed.

%once I have a list of unique states

%Figure adapted for n=4 case
iid_fig = figure;
plot(r(:,1)-1,'k') % 0 distressed
hold on
plot(r(:,9)-1,'b--') % 1 distressed
hold on
plot(r(:,10)-1,'r-.') % 2 distressed
hold on
plot(r(:,12)-1,'g.-') % 3 distressed
hold on
plot(r(:,16)-1,'c.-.') % 4 distressed
%hold on


%axis([0,T,-.04,0.13])

title(sprintf('No network, n=%d',n))
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))

legend('0 distressed','1 distressed','2 distressed',...
    '3 distressed','4 distressed',...
    'Location','northeast','Orientation','vertical')



%%
%{
cd C:\Users\Oren\Documents\MATLAB\Network
print(iid_fig,'-djpeg',...
    sprintf('iid-eta%.2f-lambda%.2f',eta,lambda))
%}








