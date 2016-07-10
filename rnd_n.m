clear; clc



%define parameters



%Total dividend for each state


n = 7; % # nodes
T = 30; % # periods forward running the risk free rate
m=12; %number of edges, maximum n(n-1)
lambda= .6;
eta=.09;

delta = .95; %time prefrence
gamma = 2; %RRA parameter of CRRA

x(1) = 1; %x_lowerbar, value of idiosyncratic component at busat state
x(2) = 2; %x_bar, value of idiosyncratic component at boom state

%Generate a random adjacency matrix of size n nodes with m edges
%NOT SYMMETRIC: decide on directionality -- incmoing links are a row
%Delta = rndgraph(n,m);
%Delta = [0,1,1,1;1,0,0,0;1,0,0,0;1,0,0,0]; %star network, two way eff
%Delta = [0,1,1,1;zeros(3,4)]'; %star network, one way eff
%Delta = zeros(n);
%Delta = [0,1,0,0;0,0,1,0;0,0,0,1;1,0,0,0]; %ring with direction
%Delta = [0,1,0,1;1,0,1,0;0,1,0,1;1,0,1,0]; %ring without direction
%Delta = ring_net(n);
%Delta = zeros(n);
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

    Lambda(:,k) = lambda*Delta*H(:,k)+lambda;
end
%%
%{
at state k, what is the prob of moving to state l?
in state l we have H(:,l)
firm i is either doing 00 01 10 11
each of those has a probability at state k
%}

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

figure;
plot(r)







%%
%{
%Demand and Supply curves for insurance as a function r_2
fig = figure;
hax=axes;
plot(r(:,1)-1,'k') % 0 periphery firm is distressed; center not
hold on
plot(r(:,2)-1,'b') % 1 periphery firm is distressed; center not
hold on
plot(r(:,4)-1,'r') % 2 periphery firm is distressed; center not
hold on
plot(r(:,8)-1,'g') % 3 periphery firm is distressed; center not
hold on
plot(r(:,9)-1,'k--') % 0 periphery firm is distressed; center yes
hold on
plot(r(:,10)-1,'b--') % 1 periphery firm is distressed; center yes
hold on
plot(r(:,12)-1,'r--') % 2 periphery firm is distressed; center yes
hold on
plot(r(:,16)-1,'g--') % 3 periphery firm is distressed; center yes

title(sprintf('Star network, n=%d',n))
xlabel(sprintf('Maturity (t)'))
ylabel(sprintf('Yield (R_f)'))

legend('No distress','1 periphery distressed',...
    '2 periphery distressed','3 periphery distressed',...
    'Only center distressed','Center & 1 periphery',...
    'Center & 2 periphery','Center & 3 periphery',...
    'Location','southeast','Orientation','vertical')
%}    
%%
%{
print(fig,'-djpeg',...
    'C:\Users\7User\Documents\MATLAB\Network\figure')
%}








