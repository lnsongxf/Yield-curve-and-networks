clear; clc

n = 5;

num_states = 2^n;

%Define states of the world (W)
W = zeros(n,2^n);
for j=1:2^n
    str = dec2bin(j-1,n);
    x = str-'0';
    W(:,j) = x';    
end

%Find pairs of states which are identical up to rotation
state_pairs = []; %Matrix of state-pairs
count=0;
for i=1:2^n
    clc
    count = count+1;
    sprintf('%d',count)
    for j=1:2^n
        if i < j
            for k=1:n-1
                if W(:,j)== circshift(W(:,i),k)
                    state_pairs = [state_pairs; i,j];                    
                end
            end
        end
    end
end
%%

row2remove = [];
for i=1:size(state_pairs,1)
    [row,col] = find(state_pairs(i,2) == state_pairs(:,1));
    row2remove = [row2remove; row];
end
row2remove = row2remove';
state_pairs(row2remove,:)=[];
unique_states = unique(state_pairs(:,1));
unique_states = [1; unique_states; 2^n] %Unique states of the world
%record(t) = size(M_,1)
%%
%{
A = zeros(2^n);

for i = 1:size(M,1)
    A(M(i,1),M(i,2)) = 1;
end
A = A+A'

%}








