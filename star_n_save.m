clear; clc



%define parameters
n = 13; % # nodes

center0 = [];
center1 = [];

for i = 1:n+1
    for k = 1:2^(n-1)
        if mod(k,100) == 0
            clc
            sprintf('%d/%d; %d/%d',i,n+1,k,2^n)    
        end

        aux = dec2bin(k-1,n); %generate binary number k-1
        digit = aux-'0'; %arrange each digit in a row vector
        if sum(digit) == i-1
            center0 = [center0; 2^(i-1), k, i-1];
        end

        aux = dec2bin(k-1+2^(n-1),n); %generate binary number k-1
        digit = aux-'0'; %arrange each digit in a row vector
        if sum(digit) == i && digit(1) == 1
            center1 = [center1; 2^(n-1)+2^(i-1), k+2^(n-1), i];
        end
    end
end


state_pairs = [center0; 2^(n-1), NaN(1), n-1; 2^(n-1)+1, NaN(1), 1; center1];
ind = find(state_pairs(:,1) == state_pairs(:,2));
state_pairs(ind,:) = [];


unique_states = [1; unique(state_pairs(:,1)); 2^n];

state_num_dist = sortrows(fliplr(unique(state_pairs(:,[1,3]),'rows')));
state_num_dist = sortrows([state_num_dist; 0, 1; n, 2^n]);


state_pairs = state_pairs(:,1:2);

% Count the number of unique states with x distressed
unique_state_x_dist = NaN(n+1,1);
for k=1:n+1
    kk=k-1;
    unique_state_x_dist(k) = sum(state_num_dist(:,1)==kk);
end

clearvars -except n state_pairs unique_states unique_state_x_dist state_num_dist

save(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\star_%d',n))