clear;

[sigma_x,sigma_y,sigma_z,I,Ix,Iy,Iz,ST0,ST1,KIx,KIy,KIz] = MultiPauli(3);

rho_in = kron(ST0,kron(ST0,ST1));
N = 600;
step = 2*pi/N;
X = 0:step:(2*pi);

%% Case 1: Symmetric XX+YY %%

ii = 1;
for t = 0:step:(2*pi)

H12 = 2*(KIx{1}*KIx{2}+KIy{1}*KIy{2});
U12 = expm(-i*H12*t);
H23 = 2*(KIx{2}*KIx{3}+KIy{2}*KIy{3});
U23 = expm(-i*H23*t);
H13 = 2*(KIx{1}*KIx{3}+KIy{1}*KIy{3});
U13 = expm(-i*H13*t);

rho_fin = (U12*U23*U13*U13*U23*U12)*rho_in*(U12*U23*U13*U13*U23*U12)';

A(ii) = abs(trace(rho_fin*kron(ST1,kron(ST0,ST0))));
B(ii) = abs(trace(rho_fin*kron(ST0,kron(ST1,ST0))));
C(ii) = abs(trace(rho_fin*kron(ST0,kron(ST0,ST1))));
D(ii) = A(ii)+B(ii)+C(ii);

ii = ii+1;

end

A_exp = [0.029 0.065 0.391 0.428 0.399 0.187 0.037 0.014 0.008 0.011...
         0.034 0.011 0.062 0.173 0.447 0.477 0.349 0.054 0.002 0.054...
         0.312 0.458 0.393 0.190 0.047 0.001 0.009 0.024 0.037 0.028...
         0.054 0.199 0.380 0.438 0.314 0.040 0.005];          
B_exp = [0.048 0.110 0.280 0.378 0.384 0.264 0.144 0.078 0.009 0.016...
         0.027 0.077 0.164 0.340 0.454 0.398 0.218 0.078 0.049 0.144...
         0.295 0.364 0.355 0.297 0.181 0.068 0.043 0.032 0.024 0.016...
         0.132 0.314 0.334 0.360 0.255 0.097 0.013];          
C_exp = [0.879 0.688 0.255 0.036 0.077 0.309 0.559 0.808 0.904 0.908...
         0.872 0.751 0.674 0.355 0.075 0.040 0.321 0.720 0.883 0.700...
         0.307 0.045 0.132 0.322 0.747 0.797 0.866 0.885 0.884 0.871...
         0.647 0.352 0.094 0.040 0.252 0.686 0.920];       
     
Y = 0:(2*pi/36):(length(A_exp)*2*pi/36-2*pi/36);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Parameters %%
subplot(2,2,1);
plot(X,real(A),'b-','Linewidth',2,'Color',[0.4667 0.7647 0.3882])
hold on
plot(X,real(B),'r--','Linewidth',2,'Color',[0.4 0.4 1])
hold on
plot(X,real(C),'g-.','Linewidth',2,'Color',[1 0.6 0.6])

plot(Y,A_exp,'b^','Marker','^','MarkerSize',6,'Linewidth',2,'Color',[0.4667 0.7647 0.3882]);
plot(Y,B_exp,'bo','Marker','o','MarkerSize',6,'Linewidth',2,'Color',[0.4 0.4 1]);
plot(Y,C_exp,'bs','Marker','s','MarkerSize',6,'Linewidth',2,'Color',[1 0.6 0.6]);

fontsize = 14;
fontsize_legend = 8;


title({'XX+YY'},'FontWeight','normal','FontSize',fontsize,...
   'FontName','Calibri',...
   'FontAngle','normal');

legend('p(3\rightarrow 1)','p(3\rightarrow 2)','p(3\rightarrow 3)');
set(legend,'EdgeColor',[1 0.8 0],'YColor',[1 0.8 0],'XColor',[1 0.8 0],...
    'LineWidth',2,...
    'FontWeight','normal','FontSize',fontsize_legend);
set(legend,'Color','none','box','off')

xlabel({'t'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');

ylabel({'p(3\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

% px=[0:pi/2:2*pi];
% py=px*0-0.05;
% text(px,py,{'0','\pi/2','\pi','3\pi/2','2\pi'},'FontWeight','normal','FontName','Calibri','FontSize',14);
% text(pi,-0.15,{'t'},'FontWeight','normal','FontName','Calibri','FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Case 2: Symmetric -(XX+YY) %%

ii = 1;
for t = 0:step:(2*pi)

H12 = -2*(KIx{1}*KIx{2}+KIy{1}*KIy{2});
U12 = expm(-i*H12*t);
H23 = -2*(KIx{2}*KIx{3}+KIy{2}*KIy{3});
U23 = expm(-i*H23*t);
H13 = -2*(KIx{1}*KIx{3}+KIy{1}*KIy{3});
U13 = expm(-i*H13*t);

rho_fin = (U12*U23*U13*U13*U23*U12)*rho_in*(U12*U23*U13*U13*U23*U12)';

A(ii) = abs(trace(rho_fin*kron(ST1,kron(ST0,ST0))));
B(ii) = abs(trace(rho_fin*kron(ST0,kron(ST1,ST0))));
C(ii) = abs(trace(rho_fin*kron(ST0,kron(ST0,ST1))));
D(ii) = A(ii)+B(ii)+C(ii);

ii = ii+1;

end

 
% A_exp = [0.011 0.059 0.229 0.207 0.157 0.012 0.082 0.156 0.109 0.003...
%          0.072 0.143 0.099 0.042 0.152 0.326 0.251 0.048 0.042 0.085...
%          0.210 0.311 0.136 0.068 0.101 0.181 0.047 0.052 0.072 0.153...
%          0.073 0.029 0.092 0.288 0.219 0.045 0.027];  
A_exp = [0.057 0.102 0.341 0.436 0.495 0.152 0.051 0.036 0.032 0.017...
         0.048 0.006 0.045 0.220 0.437 0.438 0.360 0.049 0.011 0.108...
         0.352 0.532 0.447 0.228 0.037 0.027 0.005 0.008 0.028 0.021...
         0.048 0.236 0.422 0.492 0.355 0.077 0.027];  
% B_exp = [0.899 0.697 0.321 0.192 0.437 0.610 0.652 0.718 0.776 0.885...
%          0.825 0.696 0.657 0.620 0.358 0.250 0.345 0.797 0.903 0.767...
%          0.421 0.242 0.411 0.599 0.634 0.656 0.827 0.897 0.815 0.702...
%          0.617 0.594 0.425 0.261 0.393 0.636 0.945];   
B_exp = [0.032 0.062 0.217 0.356 0.452 0.305 0.186 0.036 0.057 0.026...
         0.045 0.053 0.226 0.315 0.336 0.355 0.249 0.075 0.026 0.120...
         0.278 0.355 0.412 0.307 0.183 0.013 0.016 0.002 0.007 0.028...
         0.173 0.322 0.389 0.366 0.247 0.124 0.036];  
% C_exp = [0.022 0.074 0.267 0.379 0.453 0.278 0.112 0.033 0.006 0.024...
%          0.029 0.033 0.112 0.277 0.372 0.380 0.264 0.081 0.032 0.066...
%          0.261 0.421 0.449 0.290 0.153 0.020 0.000 0.037 0.033 0.060...
%          0.083 0.245 0.379 0.371 0.305 0.058 0.044];      
C_exp = [0.858 0.766 0.317 0.046 0.154 0.353 0.696 0.789 0.878 0.914...
         0.857 0.836 0.619 0.376 0.045 0.008 0.276 0.690 0.854 0.759...
         0.327 0.032 0.149 0.360 0.588 0.796 0.880 0.869 0.876 0.858...
         0.708 0.413 0.144 0.092 0.290 0.665 0.852];  
%      
Y = 0:(2*pi/36):(length(A_exp)*2*pi/36-2*pi/36);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Parameters %%
subplot(2,2,2);
plot(X,real(A),'b-','Linewidth',2,'Color',[0.4667 0.7647 0.3882])
hold on
plot(X,real(B),'r--','Linewidth',2,'Color',[0.4 0.4 1])
hold on
plot(X,real(C),'g-.','Linewidth',2,'Color',[1 0.6 0.6])

plot(Y,A_exp,'b^','Marker','^','MarkerSize',6,'Linewidth',2,'Color',[0.4667 0.7647 0.3882]);
plot(Y,B_exp,'bo','Marker','o','MarkerSize',6,'Linewidth',2,'Color',[0.4 0.4 1]);
plot(Y,C_exp,'bs','Marker','s','MarkerSize',6,'Linewidth',2,'Color',[1 0.6 0.6]);

title({'-(XX+YY)'},'FontWeight','normal','FontSize',fontsize,...
   'FontName','Calibri',...
   'FontAngle','normal');

% legend('p(1\rightarrow 1)','p(1\rightarrow 2)','p(1\rightarrow 3)');
% set(legend,'EdgeColor',[1 0.8 0],'YColor',[1 0.8 0],'XColor',[1 0.8 0],...
%     'LineWidth',2,...
%     'FontWeight','normal','FontSize',fontsize_legend);
%set(legend,'Color','none')

xlabel({'t'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');

ylabel({'p(3\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

% px=[0:pi/2:2*pi];
% py=px*0-0.05;
% text(px,py,{'0','\pi/2','\pi','3\pi/2','2\pi'},'FontWeight','normal','FontName','Calibri','FontSize',14);
% text(pi,-0.15,{'t'},'FontWeight','normal','FontName','Calibri','FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Case 3: AntiSymmetric XY-YX %%

ii = 1;
for t = 0:step:(2*pi)

H12 = 2*(KIx{1}*KIy{2}-KIy{1}*KIx{2});
U12 = expm(-i*H12*t);
H23 = 2*(KIx{2}*KIy{3}-KIy{2}*KIx{3});
U23 = expm(-i*H23*t);
H13 = 2*(KIx{1}*KIy{3}-KIy{1}*KIx{3});
U13 = expm(-i*H13*t);

rho_fin = (U12*U23*U13*U13*U23*U12)*rho_in*(U12*U23*U13*U13*U23*U12)';

A(ii) = abs(trace(rho_fin*kron(ST1,kron(ST0,ST0))));
B(ii) = abs(trace(rho_fin*kron(ST0,kron(ST1,ST0))));
C(ii) = abs(trace(rho_fin*kron(ST0,kron(ST0,ST1))));
D(ii) = A(ii)+B(ii)+C(ii);

ii = ii+1;

end

    
  
A_exp = [0.042 0.065 0.085 0.108 0.027 0.052 0.013 0.016 0.007 0.020...
         0.000 0.015 0.011 0.032 0.023 0.125 0.102 0.038 0.026 0.120...
         0.465 0.778 0.684 0.517 0.147 0.014 0.008 0.030 0.015 0.028...
         0.138 0.447 0.738 0.856 0.468 0.068 0.036];  
B_exp = [0.055 0.163 0.470 0.726 0.840 0.526 0.194 0.026 0.020 0.044...
         0.020 0.040 0.162 0.448 0.842 0.787 0.467 0.163 0.001 0.089...
         0.135 0.065 0.004 0.025 0.038 0.019 0.039 0.019 0.005 0.044...
         0.083 0.075 0.039 0.044 0.062 0.085 0.068];  
C_exp = [0.884 0.636 0.328 0.032 0.126 0.346 0.684 0.846 0.878 0.928...
         0.902 0.846 0.629 0.399 0.067 0.056 0.274 0.700 0.895 0.749...
         0.318 0.054 0.132 0.319 0.647 0.819 0.889 0.879 0.860 0.792...
         0.644 0.329 0.087 0.050 0.241 0.749 0.929];
% %      
Y = 0:(2*pi/36):(length(A_exp)*2*pi/36-2*pi/36);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Parameters %%
subplot(2,2,3);
plot(X,real(A),'b-','Linewidth',2,'Color',[0.4667 0.7647 0.3882])
hold on
plot(X,real(B),'r--','Linewidth',2,'Color',[0.4 0.4 1])
hold on
plot(X,real(C),'g-.','Linewidth',2,'Color',[1 0.6 0.6])

plot(Y,A_exp,'b^','Marker','^','MarkerSize',6,'Linewidth',2,'Color',[0.4667 0.7647 0.3882]);
plot(Y,B_exp,'bo','Marker','o','MarkerSize',6,'Linewidth',2,'Color',[0.4 0.4 1]);
plot(Y,C_exp,'bs','Marker','s','MarkerSize',6,'Linewidth',2,'Color',[1 0.6 0.6]);

title({'XY-YX'},'FontWeight','normal','FontSize',fontsize,...
   'FontName','Calibri',...
   'FontAngle','normal');

% legend('p(1\rightarrow 1)','p(1\rightarrow 2)','p(1\rightarrow 3)');
% set(legend,'EdgeColor',[1 0.8 0],'YColor',[1 0.8 0],'XColor',[1 0.8 0],...
%     'LineWidth',2,...
%     'FontWeight','normal','FontSize',fontsize_legend);
%set(legend,'Color','none')

xlabel({'t'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');

ylabel({'p(3\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

% px=[0:pi/2:2*pi];
% py=px*0-0.05;
% text(px,py,{'0','\pi/2','\pi','3\pi/2','2\pi'},'FontWeight','normal','FontName','Calibri','FontSize',14);
% text(pi,-0.15,{'t'},'FontWeight','normal','FontName','Calibri','FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Case 4: AntiSymmetric -(XY-YX) %%

ii = 1;
for t = 0:step:(2*pi)

H12 = -2*(KIx{1}*KIy{2}-KIy{1}*KIx{2});
U12 = expm(-i*H12*t);
H23 = -2*(KIx{2}*KIy{3}-KIy{2}*KIx{3});
U23 = expm(-i*H23*t);
H13 = -2*(KIx{1}*KIy{3}-KIy{1}*KIx{3});
U13 = expm(-i*H13*t);

rho_fin = (U12*U23*U13*U13*U23*U12)*rho_in*(U12*U23*U13*U13*U23*U12)';

A(ii) = abs(trace(rho_fin*kron(ST1,kron(ST0,ST0))));
B(ii) = abs(trace(rho_fin*kron(ST0,kron(ST1,ST0))));
C(ii) = abs(trace(rho_fin*kron(ST0,kron(ST0,ST1))));
D(ii) = A(ii)+B(ii)+C(ii);

ii = ii+1;

end

A_exp = [0.002 0.098 0.446 0.732 0.668 0.494 0.173 0.016 0.006 0.021...
         0.055 0.012 0.165 0.456 0.755 0.794 0.496 0.096 0.011 0.055...
         0.106 0.080 0.054 0.024 0.016 0.051 0.011 0.015 0.039 0.016...
         0.016 0.016 0.048 0.105 0.059 0.030 0.003];
B_exp = [0.029 0.072 0.069 0.055 0.036 0.069 0.081 0.027 0.007 0.041...
         0.038 0.017 0.041 0.054 0.027 0.020 0.145 0.045 0.015 0.102...
         0.417 0.708 0.767 0.455 0.189 0.064 0.026 0.005 0.034 0.034...
         0.183 0.382 0.763 0.750 0.477 0.155 0.027];
C_exp = [0.889 0.724 0.325 0.005 0.068 0.303 0.656 0.763 0.849 0.878...
         0.858 0.846 0.657 0.327 0.125 0.018 0.233 0.737 0.901 0.706...
         0.238 0.061 0.070 0.327 0.717 0.824 0.872 0.936 0.928 0.852...
         0.667 0.357 0.119 0.054 0.328 0.567 0.872]; 
% %      
Y = 0:(2*pi/36):(length(A_exp)*2*pi/36-2*pi/36);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Parameters %%
subplot(2,2,4);
plot(X,real(A),'b-','Linewidth',2,'Color',[0.4667 0.7647 0.3882])
hold on
plot(X,real(B),'r--','Linewidth',2,'Color',[0.4 0.4 1])
hold on
plot(X,real(C),'g-.','Linewidth',2,'Color',[1 0.6 0.6])

plot(Y,A_exp,'b^','Marker','^','MarkerSize',6,'Linewidth',2,'Color',[0.4667 0.7647 0.3882]);
plot(Y,B_exp,'bo','Marker','o','MarkerSize',6,'Linewidth',2,'Color',[0.4 0.4 1]);
plot(Y,C_exp,'bs','Marker','s','MarkerSize',6,'Linewidth',2,'Color',[1 0.6 0.6]);

title({'-(XY-YX)'},'FontWeight','normal','FontSize',fontsize,...
   'FontName','Calibri',...
   'FontAngle','normal');

% legend('p(1\rightarrow 1)','p(1\rightarrow 2)','p(1\rightarrow 3)');
% set(legend,'EdgeColor',[1 0.8 0],'YColor',[1 0.8 0],'XColor',[1 0.8 0],...
%     'LineWidth',2,...
%     'FontWeight','normal','FontSize',fontsize_legend);
%set(legend,'Color','none')

xlabel({'t'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');

ylabel({'p(3\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

% px=[0:pi/2:2*pi];
% py=px*0-0.05;
% text(px,py,{'0','\pi/2','\pi','3\pi/2','2\pi'},'FontWeight','normal','FontName','Calibri','FontSize',14);
% text(pi,-0.15,{'t'},'FontWeight','normal','FontName','Calibri','FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%