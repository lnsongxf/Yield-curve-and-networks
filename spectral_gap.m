clear; clc

range = 300;
precision = 1000;
sg_record = NaN(range);
sg_lambda_axis = NaN(range,1);
sg_eta_axis = NaN(range,1);

for sg1 = 1:range
    for sg2 = 1:range
    clc
    sprintf('sg1: %d/%d; sg2: %d/%d',sg1,range,sg2,range)    
        

        

%define parameters
n = 6; % # nodes
T = 30; % # periods forward running the risk free rate
%m=.1;
%mm = .2;
lambda= sg1/precision; %internsity of distress
eta=sg2/precision; %probability of recovery
sg_lambda_axis(sg1) = lambda;
sg_eta_axis(sg2) = eta;

delta = .95; %time prefrence
gamma = 2; %RRA parameter of CRRA

%Dividend at each state
x(1) = 1; %x_lowerbar, value of idiosyncratic component at busat state
x(2) = 2; %x_bar, value of idiosyncratic component at boom state

Delta = ring_net(n); %

H = zeros(n,2^n); %States of the world (H)
Lambda = zeros(n,2^n); %collect all lambda_i for each node at each state
Eta = zeros(n,2^n); %collect all lambda_i for each node at each state
Pi = zeros(2^n,2^n); %markov matrix
D = NaN(2^n,1); %total dividend at each state

str = cell(2^n,1);

for k=1:2^n %state of the world
    str{k} = dec2bin(k-1,n); %generate binary number k-1
    digit = str{k}-'0'; %arrange each digit in a row vector
    H(:,k) = digit'; %define this vector as state of the world k
    D(k) = x(1)*sum(H(:,k))+x(2)*(n-sum(H(:,k))); %total dividend at k

    Lambda(:,k) = lambda*Delta*H(:,k);
    Eta(:,k) = eta*Delta*H(:,k);    
    %Lambda(:,k) = lambda*Delta*H(:,k)+lambda;
end



%Find pairs of states which are identical up to rotation
state_pairs = []; %Matrix of state-pairs
count=0;
for i=1:2^n
    %clc
    count = count+1;
    %sprintf('%d',count)
    for j=1:2^n
        if i < j
            for k=1:n-1
                if H(:,j)== circshift(H(:,i),k)
                    state_pairs = [state_pairs; i,j];
                elseif H(:,j) == circshift(fliplr(H(:,i)')',k) %mirror image
                    state_pairs = [state_pairs; i,j];
                end
            end
        end
    end
end


row2remove = [];
for i=1:size(state_pairs,1)
    [row,col] = find(state_pairs(i,2) == state_pairs(:,1));
    row2remove = [row2remove; row];
end
row2remove = row2remove';
state_pairs(row2remove,:)=[];
unique_states = unique(state_pairs(:,1));
unique_states = [1; unique_states; 2^n]; %Unique states of the world
%%
%arrange unique states by number of distressed
aux = NaN(size(unique_states,1),2);
for k = 1:size(unique_states,1)
    aux(k,:) = [sum(H(:,unique_states(k))),unique_states(k)];
end

[aux_,index]=sort(aux(:,1));

%use the column indices from sort() to sort all columns of aux.
state_num_dist = aux(index,:);
unique_states = state_num_dist(:,2);

% Count the number of unique states with x distressed
unique_state_x_dist = NaN(n+1,1);
for k=1:n+1
    kk=k-1;
    unique_state_x_dist(k) = sum(state_num_dist(:,1)==kk);    
end


%%
%Constructing the transition matrix

for k_origin = 1:2^n
    for k_dest = 1:2^n

        aux = NaN(n,1);
        for i = 1:n
aux(i) = (1-Eta(i,k_origin))*H(i,k_origin)*H(i,k_dest)+Eta(i,k_origin)*H(i,k_origin)*(1-H(i,k_dest))+...
         Lambda(i,k_origin)*(1-H(i,k_origin))*H(i,k_dest)+(1-Lambda(i,k_origin))*(1-H(i,k_origin))*(1-H(i,k_dest));
        end
        Pi(k_origin,k_dest) = prod(aux);
        
    end
end


MU = @(u) u.^(-gamma);
A = Pi.*(MU(D)*(1./MU(D))')*delta;
%sort(eig(A));
aux = eig(A);
aux1 = find(abs(aux-.95) >=.01)';
sg_record(sg1,sg2) = .95-max(aux(aux1));
%MU = @(u) exp(-delta)*u.^(-gamma);
%A = Pi.*(MU(D)*(1./MU(D))');



    end
end
%%
figure;
surf(sg_lambda_axis,sg_eta_axis,real(sg_record))




