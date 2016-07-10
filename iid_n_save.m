clear; clc

%define parameters
n = 8; % # nodes
    clc
    disp(sprintf('%d/%d',n,20))

state_num = [];
for k=1:2^n %state of the world
    str = dec2bin(k-1,n); %generate binary number k-1
    digit = str-'0'; %arrange each digit in a row vector
    state_num = [state_num; k, sum(digit)];
end

state_num = sortrows(state_num,2);
state_pairs = [];
%%
for i = 2:n    
    right_entries = state_num(state_num(:,2) == i-1,1);
    state_pairs = [state_pairs; 2^(i-1)*ones(size(right_entries(2:end),1),1), right_entries(2:end), state_num(state_num(:,1)==2^(i-1),2)*ones(size(right_entries(2:end),1),1)];
end

unique_states = unique(state_pairs(:,1));
unique_states = [1; unique_states; 2^n]; % full list of unique states


%arrange unique states by number of distressed
aux = NaN(size(unique_states,1),2);
for k = 1:size(unique_states,1)
    aux(k,:) = [state_num(state_num(:,1) == unique_states(k),2), unique_states(k)];
end


[aux_,index]=sort(aux(:,1)); %sort rows

%use the column indices from sort() to sort all columns of aux.
state_num_dist = aux(index,:)

%unique_states = state_num_dist(:,2)

% Count the number of unique states with x distressed
unique_state_x_dist = NaN(n+1,1);
for k=1:n+1
    kk=k-1;
    unique_state_x_dist(k) = sum(state_num_dist(:,1)==kk);
end


clearvars -except n state_pairs unique_states unique_state_x_dist state_num_dist
%clearvars ans aux aux_ aux_H_index aux0 aux1 aux0_H_index aux1_H_index col count digit i j k kk row row2remove x dist_num H_ index unique_states1 H unique_states flag...
%    state_pairs01 state_pairs state_pairs_ state_pairs__ Eta Lambda

save(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\iid_%d',n))

