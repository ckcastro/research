function plot_appearance(gcf,filename)
% usage plot_appearance(gcf)
%      Input -     gcf : figure handle
%              filename: string of the filename
     
% Tighthen white space from figure.
% Save pdf files. 

set(0,'defaultlinelinewidth',2) % line thickness everywhere
set(0,'defaultaxesfontsize',12)

set(0,'DefaultTextFontSize',14); % Global text fontsize

% global legend fontsize
%set(0,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual')


%legendfontsize = 12;

set(gca, 'FontSize', 14); % Axis thicks fontsize


set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);

filename = strcat(pwd,'/figures/',filename);
print('-dpdf',filename);
%print('-depsc',filename)
%print('-djpeg',filename)

%savefig(filename)


end

