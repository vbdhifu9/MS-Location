clear all;clc;close all
%% 
nn=301;
mm=301;
dx=5;
dt=0.0003;
nnn=[1:41];

load VX;
rec2=SSh(:,nnn);
k0=size(rec2,1);
s0=size(rec2,2);

%% 
load VSPTTP41;
TTP = reshape(TTP, [], size(TTP, 3)); 
TTP=TTP(:,nnn);

%% 
sta=120;
lta=2*sta;
SLTA=zeros(size(rec2));
for rr=1:s0
    for i=lta+1:k0-sta
        sum_sta=sum(abs(rec2(i:i+sta,rr)))/sta;
        sum_lta=sum(abs(rec2(i-lta:i,rr)))/lta+abs(min(rec2(:)))*0.00001;
        SLTA(i,rr)=sum_sta/sum_lta;
    end
        SLTA(:,rr)=SLTA(:,rr)/max(SLTA(:,rr));
end
ind=(isnan(SLTA));
SLTA(ind)=0;

rec2=SLTA;
nnt=size(rec2,1);

%% 
t0=450;
recm1=zeros(nn*mm,1);
%%
tic
TTP=reshape(TTP,nn,mm,41);
recm1=ones(nn,mm);
for ii=1:s0% 
     for i=1:nn%
        for j=1:mm  
            nt=round(TTP(i,j,ii)/dt)+t0;%+ttt      
            nt=min(max(nt,1),nnt);                       
%             recm1(i,j)=recm1(i,j)+rec2(nt,ii);
            recm1(i,j)=recm1(i,j)+rec2(nt,ii);
        end
    end
end


%% 
[Zd,Td]=meshgrid((0:300)*dx,(0:300)*dx);

figure
imagesc(recm1/(max(recm1(:))));
set(gca,'ydir','reverse','fontsize',16,'FontName', 'Times New Roman');   
xlabel('X (m)','fontsize',16);ylabel('Z (m)','fontsize',16,'FontName', 'Times New Roman');
set(gca, 'xtick', [1:60:301], 'xticklabel', [0 300 600 900 1200 1500]);
set(gca, 'ytick', [1:60:301], 'yticklabel', [0 300 600 900 1200 1500]);
axis([1 301 1 301]);
hc=colorbar;
set(hc,'ytick',[0 0.5 1])
set(gcf,'pos',[100 100 500 350]);
caxis([0 1]);
colormap(jet)
text(10, 16, '(b)', 'FontSize', 16,'Color', 'w','FontName', 'Times New Roman');
[N,rec_n]=size(SSh);
dt=0.0003;


figure
wigb(SSh)
set(gca,'FontSize',16)
xlabel('Trace number','FontSize',16)
ylabel('T (s)','FontSize',16)
set(gca,'ytick',1:round(N/5):N,'yticklabel',(0:round(N/5):N-1)*dt);
title('','FontSize',16)

[maxValue, maxIndex] = max(recm1(:));
[xSource, ySource] = ind2sub(size(recm1), maxIndex);
sourceImagingValue = recm1(xSource, ySource);
x = (1:nn) * dx;
y = (1:mm) * dx;

x9=x';
y9=y';
maxImagingValue = max(recm1(:));
minImagingValue = min(recm1(:));
x1=(recm1(xSource, :) - minImagingValue) / (maxImagingValue - minImagingValue);
z=(recm1(:, ySource) - minImagingValue) / (maxImagingValue - minImagingValue);
x2=x1';

figure;
plot(x, x1, 'LineWidth', 2);
xlabel('X (m)','fontsize',16);
ylabel('Normalized amplitude','fontsize',16);
title('X Profile of Normalized Imaging Value');
xlim([0 1500]);
set(gca,'fontsize',16);   

figure;
plot(y, z, 'LineWidth', 2);
xlabel('Z (m)','fontsize',16);
ylabel('Normalized amplitude','fontsize',16);
title('Z Profile of Normalized Imaging Value');
xlim([0 1500]);
set(gca,'fontsize',16);  









