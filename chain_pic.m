function [ pic ] = k_reg_pic( n, s_number)

%parameters defined by function
%n = 8; % # nodes
%s_number = 13;

step = (n-1)/n*pi;
network = [linspace(-step,step,n)', zeros(n,1)];


%state of the world s_number
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

%line
plot(network([1,end],1),network([1,end],2),'-k','MarkerSize',ms,'LineWidth',lw);
hold on
%nodes
plot(network(:,1),network(:,2),'ob','MarkerSize',ms,'LineWidth',lw);
hold on
plot(network(ind,1),network(ind,2),'*r','MarkerSize',ms,'LineWidth',lw);

base = step;
base_addition = 1;
border = base+base_addition;
axis([-border,border,-.1,.1]);
axis off;



end
