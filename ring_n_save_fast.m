%{
THE PURPOSE OF THIS FILE IS TO GENERATE A LIST OF UNIQUE STATES OF THE
WORLD IN A BINARY CIRCLE. 

IN AN N-NODES CIRCLE, A STATE IS SERIES OF 0'S AND 1'S LIKE SO
101100
THE FIRST 1 AND THE LAST 0 ARE ADJACENT
THE LATTER IS SYMMETRIC UP TO ROATION TO
010110
001011
100101
110010
ETC
AND UP TO MIRROR IMAGE AND ROTATION WITH
001101
ETC 
%}

clear; clc

% There are n nodes in the circle, each is either 1 or 0
n = 4; % # nodes

% Delta is the adjacency matrix
Delta = ring_net(n); %

H = NaN(n,2^n); %States of the world H(:,i) for state i
str = cell(2^n,1);

for k=1:2^n %state of the world
    if mod(k,100) == 0
    clc
    sprintf('%d/%d',k,2^n)    
    end
    
    str{k} = dec2bin(k-1,n); %generate binary number k-1
    digit = str{k}-'0'; %arrange each digit in a row vector
    H(:,k) = logical(digit'); %define this vector as state of the world k
end


%{
*************************************************************************
                        COMPLEMENTS
*************************************************************************
Compute for each state it's complement. E.g. 011001 <--> 100110
The way to do it: run from 1 to ceil((n+1)/2) since we look at distress
number from 0 to n/2 or immediately below n/2 if n is odd

Note! On col 1 we have all states with # of dist firms up to n/2. This will
be useful for the next stage

store it at state_pairs01
%}

state_pairs01 = []; % 2^n-1 by 2 matrix of state-pairs
counter = 0;
state_pairs01 = NaN(2^(n-1),4);
for dn = 1:floor(n/2)+1
    % dn is a dummy var
    dist_num = dn-1; % dist_num is # distressed
    aux0_H_index = find(sum(H) == dist_num); % find indices of dist_num
    % Find for each state it's equivalent
    if dist_num == n/2
        for i=1:size(aux0_H_index,2)/2
            state_pairs01(i+counter,:) = [aux0_H_index(i),2^n+1-aux0_H_index(i),dist_num,n-dist_num];
        end
    else
        for i=1:size(aux0_H_index,2)
            state_pairs01(i+counter,:) = [aux0_H_index(i),2^n+1-aux0_H_index(i),dist_num,n-dist_num];
        end        
    end
    
    counter = counter+size(aux0_H_index,2);
end

state_pairs01 = sortrows(state_pairs01,1);

%{
*************************************************************************
            Second symmetry: rotation and mirror image
*************************************************************************
%}


state_pairs_ = []; %Matrix of state-pairs
for dist_num = 1:floor(n/2) % can use only until n/2, the rest is by 0-1

    state_pairs = []; %Matrix of state-pairs
    
% if n is odd or dist_num not n/2 can use H_
    aux_H_index = state_pairs01(state_pairs01(:,3) == dist_num,1);
    if dist_num == n/2
        aux_H_index = [aux_H_index; state_pairs01(state_pairs01(:,4) == dist_num,2)];
    end
    aux = H(:,aux_H_index);

%Find pairs of states which are identical up to rotation and mirror
    for i=1:size(aux,2)
        clc
        sprintf('%d/%d, dist_num=%d',i,size(aux,2),dist_num)
        flag = 0;
        for j=1:size(aux,2)
            if i < j
                counter = 0;
                while counter < n
                    if isequal(aux(:,i), circshift(aux(:,j),[counter,0])) == 1 %simple rotation
                        state_pairs = [state_pairs; ...
                            aux_H_index(i),aux_H_index(j)];
                        flag = 1;
                        counter = n;
                    elseif isequal(aux(:,i), circshift(fliplr(aux(:,j)')',[counter,0])) == 1 %mirror image
                        state_pairs = [state_pairs; ...
                            aux_H_index(i),aux_H_index(j)];
                        flag = 1;
                        counter = n;
                    end
                    counter = counter+1;
                end
            end
        end

        if flag == 0 % the state is unique
            state_pairs = [state_pairs; aux_H_index(i), NaN(1)];
        end
    end
        
        

    
    state_pairs_ = [state_pairs_; unique(state_pairs,'rows')]; %Matrix of state-pairs


end

% remove duplicates
row2remove = [];
for i=1:size(state_pairs_,1)
    [row,col] = find(state_pairs_(:,1) == state_pairs_(i,2));
    row2remove = [row2remove; row];
end

row2remove = unique(row2remove)';
state_pairs_(row2remove,:)=[];
state_pairs_ = sortrows(state_pairs_);



%{
So far we have eliminated rows that are duplicates. But some states on the
left column of state_pairs_ may have an NaN on the right -- the following
part deals with those.
%}

%find indices of unique states that have no match
ind = find(isnan(state_pairs_(:,2))); %ind in state_pairs_
ind_ = state_pairs_(isnan(state_pairs_(:,2)),1); %num of state

row2remove = [];
for i = 1:size(ind,1)
%if there is more than one entry of that unique state then remove the row
%with the NaN on the right
    if sum(state_pairs_(:,1) == ind_(i)) > 1 
        row2remove = [row2remove; ind(i)];
    end    
end

state_pairs_(row2remove,:)=[];

%{
we are now left with a list of state_pairs where on the left col we
have unique set of states and on the right col states that correspond to
them in H. But this is done only for half of the matrix -- now we find the
list of complements. Need to be careful of states with dist=n/2
%}

% find the translation of the latter list into their complements
state_pairs__ = [];
for i = 1:size(state_pairs_,1)
	if sum(H(:,state_pairs_(i,1))) ~= n/2
        state_pairs__(i,1) = ...
        state_pairs01(state_pairs01(:,1) == state_pairs_(i,1),2);

    
%some states are unique (have NaN on the right), so then when trying to
%find a complement there would not be one. could just add NaN(1,2) to the
%list of complements - both would work
        if isnan(state_pairs_(i,2)) == 0
            state_pairs__(i,2) = ...
        state_pairs01(state_pairs01(:,1) == state_pairs_(i,2),2);
        else
            state_pairs__(i,2) = NaN(1);
        end
    end
end

state_pairs__(state_pairs__(:,1)==0,:) = [];
state_pairs = [state_pairs_; state_pairs__]; % full list of state pairs
unique_states = unique(state_pairs(:,1));
unique_states = [1; unique_states; 2^n]; % full list of unique states

%arrange unique states by number of distressed
aux = NaN(size(unique_states,1),2);
for k = 1:size(unique_states,1)
    aux(k,:) = [sum(H(:,unique_states(k))),unique_states(k)];
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


clearvars -except n state_pairs unique_states state_pairs01 unique_state_x_dist state_num_dist
%clearvars ans aux aux_ aux_H_index aux0 aux1 aux0_H_index aux1_H_index col count digit i j k kk row row2remove x dist_num H_ index unique_states1 H unique_states flag...
%    state_pairs01 state_pairs state_pairs_ state_pairs__ Eta Lambda

save(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\ring_%d',n))


