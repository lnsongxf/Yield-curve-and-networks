function [ pic ] = k_reg_pic( n, s_number )
%%
%parameters defined by function
%n = 8; % # nodes
%reg = 4;
%s_number = 129;


center = [0,0]; %center of ring
r = [4,4]';
aux = [];
for i = 1:n-1;
    aux = [aux; r(mod(i-1,2)+1)];
end
r = aux;

theta = linspace(0,(n-2)/(n-1)*2*pi,n-1)'+pi/2;
network = [r.*cos(theta),r.*sin(theta); center];


%state of the world s

str_ = dec2bin(s_number-1,n);
digit = str_-'0';
s_logicals = fliplr(logical(digit))';
ind = find(s_logicals);


%set(h, 'Position', [x_pos y_pos width_fig height_fig])
ms = 7;
lw = 2;
figure;
set(gcf, 'Visible', 'off');
pic = subplot(1,1,1);

hold on
plot(network(:,1),network(:,2),'ob','MarkerSize',ms,'LineWidth',lw);
for i = 1:n
    hold on
    plot([network(n,1),network(mod(i,n)+1,1)],[network(n,2),network(mod(i,n)+1,2)],'-k')
end

hold on
plot(network(ind,1),network(ind,2),'*r','MarkerSize',ms,'LineWidth',lw);


base = r(1);
base_addition = 1;
border = base+base_addition;
axis([-border,border,-border,border]);
axis off;


%%
end
