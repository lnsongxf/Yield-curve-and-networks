clear; clc

for n = 4:8
    filename='temp.mat';    
    load(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\ring_%d',n),...
        'state_pairs', 'unique_states', 'n')
    Ring(1,1:3) = [2, 3, 4];
    
	Ring(1,n) = size(unique_states,1);
    save(filename,'Ring','-append');
end


for n = 9:16
    load(sprintf('C:\\Users\\Oren\\Documents\\MATLAB\\Network\\data\\ring_%d',n),...
        'state_pairs', 'unique_states', 'n')
    filename='temp.mat';
	Ring(2,n-8) = size(unique_states,1);
    save(filename,'Ring','-append');
end

N = [1:8; 9:16];
Star = [2:2:16; 18:2:32];
Max = NaN(2,8);
for i = 1:8
    Max(:,i) = [2^i; 2^(i+8)];
end

data = [N(1,:); Star(1,:); Ring(1,:); Max(1,:);...
        N(2,:); Star(2,:); Ring(2,:); Max(2,:)];

Rows = strread(num2str(1:8),'%s')';

T1 = array2table(data())
T2 = array2table(data)
%{
T1 = table(data(1,:)',data(2,:)',data(3,:)','RowNames',Rows)
T2 = table(data(4,:)',data(5,:)',data(6,:)','RowNames',Rows)

T1.Properties.VariableNames = {'Star','Ring','Max'}
T2.Properties.VariableNames = {'Star','Ring','Max'}
%}

filename = 'temp.csv';
csvwrite(filename,data)
%writetable(T,filename,'WriteVariableNames',false)




