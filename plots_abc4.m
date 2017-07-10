% Claudia Castro-Castro
% Fall 2015 
% Plot second moment, energy and norm relative errors from ABC4 symplectic
% integrator

clear all
clc
close all

delete *.pdf
delete *.eps
delete *.jpg
disp('======================================')
disp('    Plot results from ABC4 ')
disp('======================================')

% nonlinearity strength < beta >
betas = load('betas_file.dat')

% step size < h or tau or delta >
delta = load('delta_file.dat');
fprintf('step size = %g\n',delta(1));

% final time
zf = load('zfinal_file.dat');
fprintf('final time = %g\n',zf);

% Number of fibers
numfib = load('numfib_file.dat');
fprintf('lattice size = %g\n',numfib);

% middle fiber
middlefiber =  ( numfib + 1 )/ 2 

time = load('timefile.dat');
log10time=log10(time);

% how many betas 
ib = length(betas);

% graphs' color code < depending on how many betas >
%plotStyle = ['m','b','g','r'];


plotStyle = [      0    0.4470    0.7410;
              0.8500    0.3250    0.0980;
              0.9290    0.6940    0.1250;
              0.4940    0.1840    0.5560;
              0.4660    0.6740    0.1880;
              0.3010    0.7450    0.9330;
              0.6350    0.0780    0.1840 ];
%legendInfo = zeros(6,1);

% allocate memory
m2=zeros(length(time),ib);
Ener_rel_err=zeros(length(time),ib);
Norm_rel_err=zeros(length(time),ib);
finalNormDensity=zeros(numfib,ib);

% load files 
for j=1:length(betas)
    %m2(:, j) = load(strcat('m2_abc4_beta_',num2str(betas(j)),'.dat'));
    %par(:,j)=load(strcat('participation_abc4_beta_',num2str(beta(j)),'.dat'));
    Ener_rel_err(:,j) = load(strcat('Eerr_abc4_beta_',num2str(betas(j)),'.dat'));
    Norm_rel_err(:,j) = load(strcat('Serr_abc4_beta_',num2str(betas(j)),'.dat'));
    finalNormDensity(:,j) = load(strcat('Norm_abc4_beta_',num2str(betas(j)),'_zf.dat'));
end

numfib=length(finalNormDensity)

% Constants of motion
% Energy
figure; 
hold on;

for j=1:length(betas)
    plot(log10(time),log10(Ener_rel_err(:,j)),'Color', plotStyle(j,:),'LineWidth',1)
    legendInfo{j} = ['nonlinearity = ' num2str(betas(j))];
end
%title(['Relative error Energy  h = ' num2str(delta(1))]);
xlabel('log10 z'); ylabel('log10 Energy'); axis tight
text(0.01,0.90,'a)','Units','normalized'); 
legend(legendInfo,'Location','Best')
legend boxoff

hold off


plot_appearance(gcf,'energyrelerror_abc4')


% Norm
figure; 
hold on;
for j=1:length(betas)
    plot(log10(time),log10(Ener_rel_err(:,j)),'Color', plotStyle(j,:),'LineWidth',1)
end

%title(['Relative error Norm  h = ' num2str(delta(1))]);
%legend(['h = ' num2str(beta(j))],'Location','Best');
xlabel('log10 z'); ylabel('log10 Norm');
axis tight
text(0.01,0.90,'b)','Units','normalized'); 
legend(legendInfo,'Location','Best')
legend boxoff

hold off

plot_appearance(gcf,'normrelerror_abc4.pdf')


%% Norm density distribution at final time
ll=1:numfib;

figure; 
hold on;
for j=1:length(betas)
    plot(ll,finalNormDensity(:,j),'Color', plotStyle(j,:),'LineWidth',1)
    legend('')
end

title(['Norm density at final time = ',num2str(zf)]);
xlabel('sites'); ylabel('|\psi_n|^2');
axis tight
legend(legendInfo,'Location','Best')
legend boxoff
hold off
plot_appearance(gcf,'finalNormDensity_abc4')

% % logarithmic scale 
% figure; 
% hold on;
% for j=1:length(betas)
%     semilogy(ll,finalNormDensity(:,j),'Color', plotStyle(j,:),'LineWidth',1)
%     legend('')
% end
% 
% title(['Norm density at final time = ',num2str(zf)]);
% xlabel('sites'); ylabel('|\psi_n|^2');
% axis tight
% legend(legendInfo,'Location','Best')
% legend boxoff
% hold off
%plot_appearance(gcf,'finalNormDensitylog_abc4')



% %% m2
% 
% figure; 
% hold on;
% for j=1:length(betas)
%     %figure,
%     plot(log10(time),log10(m2(:,j)),plotStyle(j),'LineWidth',1)
%     legendInfo{j} = ['beta = ' num2str(betas(j))];
% end % loop over betas
% 
% %   Polyfit lasta beta < beta = 1 >
% cutm2 = 30;
% newtime=time(cutm2 :length(time)-1,1);
% newm2=m2(cutm2 :length(time)-1,length(betas)-1);
% 
% % fit in log domain
% polm2=polyfit(log10(newtime),log10(newm2),1);
% 
% % compute fit in linear domain
% y_hat=10.^( polm2(1)*log10(newtime) + polm2(2)); 
% %figure,
% plot(log10(newtime),log10(y_hat),'--k','LineWidth',1)  % make loglog plot
% legendInfo{j+1} = ['log10(m2) = ' num2str(polm2(1)) ' log10(t) + ' num2str(polm2(2))];
% 
% % %   Polyfit lasta beta < beta = 4.5 >
% % cutm2 =20;
% % newtime=time(cutm2 :length(time)-1,1);
% % newm2=m2(cutm2 :length(time)-1,length(betas));
% % 
% % % fit in log domain
% % polm2=polyfit(log10(newtime),log10(newm2),1);
% % 
% % % compute fit in linear domain
% % y_hat=10.^( polm2(1)*log10(newtime) + polm2(2)); 
% % %figure,
% % plot(log10(newtime),log10(y_hat),':k','LineWidth',1)  % make loglog plot
% % legendInfo{j+2} = ['log10(m2) = ' num2str(polm2(1)) ' log10(t) + ' num2str(polm2(2))];
% 
% %y=10.^( 1/3*log10(time) + polm2(2)); 
% %plot(log10(time),log10(y),'--k','LineWidth',2)
% % label2=['log10(m2) =  1/3 log10(t) + ' num2str(polm2(2))];
% %xlim([0 log10time(length(time))])
% 
% title(['Second moment ABC4 h = ' num2str(delta(1))]);
% xlabel('log10 t'); ylabel('log10 m_2');
% legend(legendInfo,'Location','Best')
% %axis([0 max(log10(time)) -50 10])
% hold off
% 
% % Make figure boundaries tight:
% ti = get(gca,'TightInset');
% set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% 
% % Now you have a tight figure on the screen but if you directly do saveas 
% % (or print), MATLAB will still add the annoying white space. To get rid of them, 
% % we need to adjust the ``paper size":
% set(gca,'units','centimeters')
% pos = get(gca,'Position');
% ti = get(gca,'TightInset');
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% 
% saveas(gcf,'m2_abc4.pdf')
% %saveas(gcf,strcat('m2_abc4.jpg'))
% %saveas(gcf,strcat('m2_abc4.eps'),'epsc')
% plot_appearance(gcf,'m2_abc4')


%% ------------------------------------------------------------------------
%     Space time contour plot
load('Normdata.dat')
fprintf(' Size normdata = '); size(Normdata)

%%

figure,
contour(Normdata(:,1:size(Normdata,2)))
%contour(Normdata)
colorbar
xlabel('propagation distance')
ylabel('lattice site')
legend boxoff
set(gca,'XLim',[0 length(time)]);% This automatically sets the
                                   % XLimMode to manual.
% Set XTick so that .
set(gca,'XTick',[0:length(time):length(time)])  % This automatically sets 
                                                  % the XTickMode to manual.    
set(gca,'XTickLabel',{'0'; num2str(zf)})
plot_appearance(gcf,'NormDensitySpaceTime')


figure,
surf(Normdata(:,1:size(Normdata,2)))


%% short time plot
figure,
contour(Normdata(:,1: ceil(size(Normdata,2)/10)))
%contour(Normdata)
colorbar
xlabel('propagation distance')
ylabel('lattice site')
legend boxoff

set(gca,'XLim',[0 ceil(size(Normdata,2)/10)]);% This automatically sets the
                                   % XLimMode to manual.
% Set XTick so that .
set(gca,'XTick',[0:ceil(size(Normdata,2)/10):ceil(size(Normdata,2)/10)])  % This automatically sets 
                                                  % the XTickMode to manual.    
set(gca,'XTickLabel',{'0'; num2str(time(ceil(size(Normdata,2)/10)))})
plot_appearance(gcf,'NormDensitySpaceShortTime')

%figure,
%surf(Normdata(:,1: ceil(size(Normdata,2)/4)))

%%
fiberZoomWindow = 50;

%NormdataZoom = Normdata(middlefiber-fiberZoomWindow:middlefiber+fiberZoomWindow,:);

zoom1=200;
zoom2=600;

NormdataZoom = Normdata(zoom1:zoom2,:);



Nc = ( length(NormdataZoom) + 1 )/ 2 

figure,
contour(NormdataZoom(:,:),200)
colorbar
xlabel('propagation distance')
ylabel('lattice site')
legend boxoff

set(gca,'YLim',[zoom1 zoom2])
%set(gca,'YTick',[zoom1:zoom1+(zoom2-zoom1)/2:zoom2])
%set(gca,'YTickLabel',[num2str(zoom1);num2str(zoom1+(zoom2-zoom1)/2);num2str(zoom1)])
%set(gca,'YTickLabel',{zoom1:zoom2}) 

set(gca,'XLim',[0 length(time)]);% This automatically sets the
                                   % XLimMode to manual.
% Set XTick so that .
set(gca,'XTick',[0:length(time):length(time)])  % This automatically sets 
                                                  % the XTickMode to manual.    
set(gca,'XTickLabel',{'0'; num2str(zf)})

plot_appearance(gcf,'NormDensitySpaceTimeZoom')

%%

NormdataZoomShortTime = Normdata(zoom1:zoom2,1:ceil(size(Normdata,2)-1)/2);
figure,
contour(NormdataZoomShortTime(:,:),200)
colorbar
xlabel('propagation distance')
ylabel('lattice site')
legend boxoff

set(gca,'YLim',[zoom1 zoom2])
%set(gca,'YTick',[zoom1:zoom1+(zoom2-zoom1)/2:zoom2])
%set(gca,'YTickLabel',[num2str(zoom1);num2str(zoom1+(zoom2-zoom1)/2);num2str(zoom1)])
%set(gca,'YTickLabel',{zoom1:zoom2}) 

set(gca,'XLim',[0 size(NormdataZoomShortTime,2)]);% This automatically sets the
                                   % XLimMode to manual.
% Set XTick so that .
set(gca,'XTick',[0:size(NormdataZoomShortTime,2):size(NormdataZoomShortTime,2)])  % This automatically sets 
                                                  % the XTickMode to manual.    
set(gca,'XTickLabel',{'0'; num2str(time(size(NormdataZoomShortTime,2)))})

plot_appearance(gcf,'NormDensitySpaceTimeZoomShortTime')


%figure,
%surf(NormdataZoom(:,:))

%% zoom and short time 
zoom1=200;
zoom2=600;

NormdataZoom = Normdata(zoom1:zoom2,ceil(size(Normdata,2)/2));

Nc = ( length(NormdataZoom) + 1 )/ 2 

figure,
contour(NormdataZoom(:,:))
colorbar
xlabel('propagation distance')
ylabel('lattice site')
legend boxoff

set(gca,'YLim',[zoom1 zoom2])
%set(gca,'YTick',[zoom1:zoom1+(zoom2-zoom1)/2:zoom2])
%set(gca,'YTickLabel',[num2str(zoom1);num2str(zoom1+(zoom2-zoom1)/2);num2str(zoom1)])
%set(gca,'YTickLabel',{zoom1:zoom2}) 

set(gca,'XLim',[0 ceil(size(Normdata,2)/2)]);% This automatically sets the
                                   % XLimMode to manual.
% Set XTick so that .
set(gca,'XTick',[0:ceil(size(Normdata,2)/2):ceil(size(Normdata,2)/2)])  % This automatically sets 
                                                  % the XTickMode to manual.    
set(gca,'XTickLabel',{'0'; num2str(time(ceil(size(Normdata,2)/2)))})

