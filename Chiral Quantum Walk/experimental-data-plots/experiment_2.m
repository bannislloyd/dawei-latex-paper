clear;

[sigma_x,sigma_y,sigma_z,I,Ix,Iy,Iz,ST0,ST1,KIx,KIy,KIz] = MultiPauli(3);

rho_in = kron(ST0,kron(ST1,ST0));
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

A_exp = [0.061 0.092 0.242 0.336 0.078 0.031 0.078 0.152 0.020 0.038...
         0.063 0.159 0.100 0.051 0.076 0.222 0.223 0.071 0.041 0.119...
         0.267 0.230 0.067 0.029 0.087 0.154 0.066 0.055 0.076 0.157...
         0.092 0.034 0.073 0.217 0.281 0.085 0.038];     
B_exp = [0.940 0.756 0.345 0.226 0.451 0.602 0.637 0.665 0.832 0.942...
         0.814 0.654 0.626 0.612 0.437 0.189 0.379 0.664 0.873 0.640...
         0.324 0.195 0.404 0.663 0.677 0.711 0.831 0.872 0.815 0.717...
         0.622 0.596 0.453 0.186 0.310 0.709 0.910];     
C_exp = [0.018 0.069 0.285 0.410 0.421 0.274 0.132 0.035 0.018 0.044...
         0.044 0.017 0.094 0.227 0.352 0.377 0.317 0.057 0.024 0.084...
         0.323 0.355 0.389 0.300 0.096 0.046 0.024 0.032 0.042 0.037...
         0.089 0.224 0.459 0.353 0.291 0.073 0.040];     
     
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

legend('p(2\rightarrow 1)','p(2\rightarrow 2)','p(2\rightarrow 3)');
set(legend,'EdgeColor',[1 0.8 0],'YColor',[1 0.8 0],'XColor',[1 0.8 0],...
    'LineWidth',2,...
    'FontWeight','normal','FontSize',fontsize_legend);
set(legend,'Color','none','box','off')

xlabel({'t'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');

ylabel({'p(2\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
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

 
A_exp = [0.011 0.059 0.229 0.207 0.157 0.012 0.082 0.156 0.109 0.003...
         0.072 0.143 0.099 0.042 0.152 0.326 0.251 0.048 0.042 0.085...
         0.210 0.311 0.136 0.068 0.101 0.181 0.047 0.052 0.072 0.153...
         0.073 0.029 0.092 0.288 0.219 0.045 0.027];     
B_exp = [0.899 0.697 0.321 0.192 0.437 0.610 0.652 0.718 0.776 0.885...
         0.825 0.696 0.657 0.620 0.358 0.250 0.345 0.797 0.903 0.767...
         0.421 0.242 0.411 0.599 0.634 0.656 0.827 0.897 0.815 0.702...
         0.617 0.594 0.425 0.261 0.393 0.636 0.945];      
C_exp = [0.022 0.074 0.267 0.379 0.453 0.278 0.112 0.033 0.006 0.024...
         0.029 0.033 0.112 0.277 0.372 0.380 0.264 0.081 0.032 0.066...
         0.261 0.421 0.449 0.290 0.153 0.020 0.000 0.037 0.033 0.060...
         0.083 0.245 0.379 0.371 0.305 0.058 0.044];      
     
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

ylabel({'p(2\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
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

    
A_exp = [0.038 0.189 0.477 0.718 0.845 0.765 0.487 0.263 0.114 0.043...
         0.069 0.282 0.508 0.695 0.827 0.744 0.458 0.134 0.029 0.088...
         0.059 0.038 0.095 0.261 0.419 0.260 0.094 0.037 0.079 0.266...
         0.369 0.300 0.067 0.002 0.054 0.055 0.034];  
B_exp = [0.903 0.651 0.384 0.156 0.039 0.069 0.290 0.544 0.782 0.936...
         0.788 0.531 0.341 0.102 0.106 0.159 0.419 0.765 0.882 0.677...
         0.304 0.101 0.041 0.077 0.292 0.518 0.773 0.895 0.792 0.509...
         0.327 0.093 0.044 0.178 0.293 0.722 0.868];  
C_exp = [0.042 0.036 0.059 0.047 0.039 0.031 0.070 0.021 0.032 0.036...
         0.017 0.021 0.030 0.044 0.037 0.074 0.061 0.071 0.001 0.116...
         0.384 0.734 0.689 0.455 0.201 0.022 0.034 0.020 0.027 0.025...
         0.178 0.441 0.682 0.768 0.558 0.085 0.031];
%      
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

ylabel({'p(2\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
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

 
A_exp = [0.049 0.041 0.083 0.007 0.085 0.267 0.365 0.262 0.139 0.012...
         0.115 0.285 0.425 0.293 0.117 0.022 0.122 0.028 0.019 0.167...
         0.434 0.718 0.804 0.736 0.556 0.226 0.062 0.011 0.070 0.279...
         0.461 0.715 0.813 0.672 0.439 0.125 0.003];   
B_exp = [0.849 0.667 0.374 0.163 0.086 0.062 0.216 0.475 0.689 0.904...
         0.809 0.406 0.234 0.124 0.034 0.085 0.307 0.668 0.926 0.787...
         0.460 0.151 0.075 0.067 0.269 0.585 0.785 0.861 0.728 0.530...
         0.292 0.067 0.077 0.081 0.429 0.648 0.848]; 
C_exp = [0.020 0.110 0.493 0.754 0.705 0.470 0.166 0.023 0.003 0.054...
         0.015 0.025 0.161 0.462 0.678 0.703 0.489 0.082 0.009 0.082...
         0.072 0.018 0.020 0.065 0.046 0.015 0.005 0.010 0.003 0.017...
         0.097 0.022 0.029 0.049 0.064 0.031 0.008]; 
%      
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

ylabel({'p(2\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

% px=[0:pi/2:2*pi];
% py=px*0-0.05;
% text(px,py,{'0','\pi/2','\pi','3\pi/2','2\pi'},'FontWeight','normal','FontName','Calibri','FontSize',14);
% text(pi,-0.15,{'t'},'FontWeight','normal','FontName','Calibri','FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%