clear; clc

%define parameters
n = 8; % # nodes

aux = [];
for k=1:2^n %state of the world
    str = dec2bin(k-1,n); %generate binary number k-1
    digit = str-'0'; %arrange each digit in a row vector
    aux = [aux; k, sum(digit)];
end

aux = sortrows(aux,2);
state_pairs = [];
%%
for i = 2:n    
    right_entries = aux(aux(:,2) == i-1,1);
    state_pairs = [state_num_dist; 2^(i-1)*ones(size(right_entries(2:end),1),1), right_entries(2:end), aux(aux(:,1)==2^(i-1),2)*ones(size(right_entries(2:end),1),1)];
end

state_pairs = [1, NaN(1)

