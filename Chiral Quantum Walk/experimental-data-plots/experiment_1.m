clear;

[sigma_x,sigma_y,sigma_z,I,Ix,Iy,Iz,ST0,ST1,KIx,KIy,KIz] = MultiPauli(3);

rho_in = kron(ST1,kron(ST0,ST0));
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



A_exp = [0.899 0.692 0.321 0.148 0.386 0.616 0.676 0.689 0.839 0.892...
         0.791 0.698 0.681 0.653 0.363 0.152 0.336 0.685 0.914...
         0.596 0.257 0.178 0.359 0.573 0.693 0.624 0.788 0.964...
         0.770 0.691 0.666 0.648 0.355 0.148 0.405 0.591 0.904];
B_exp = [0.017 0.122 0.249 0.193 0.110 0.032 0.048 0.192 0.102 0.053...
         0.099 0.163 0.110 0.025 0.048 0.286 0.238 0.096 0.046...
         0.084 0.122 0.287 0.115 0.021 0.073 0.140 0.062 0.053...
         0.167 0.224 0.122 0.014 0.033 0.270 0.162 0.170 0.022];
C_exp = [0.009 0.083 0.312 0.441 0.365 0.234 0.051 0.034 0.027 0.034...
         0.010 0.036 0.054 0.225 0.439 0.476 0.294 0.064 0.029...
         0.053 0.256 0.396 0.438 0.195 0.104 0.064 0.035 0.052...
         0.036 0.030 0.100 0.155 0.384 0.434 0.289 0.080 0.038];
     
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

legend('p(1\rightarrow 1)','p(1\rightarrow 2)','p(1\rightarrow 3)');
set(legend,'EdgeColor',[1 0.8 0],'YColor',[1 0.8 0],'XColor',[1 0.8 0],...
    'LineWidth',2,...
    'FontWeight','normal','FontSize',fontsize_legend);
set(legend,'Color','none','box','off')

xlabel({'t'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');

ylabel({'p(1\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

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


A_exp = [0.916 0.710 0.350 0.093 0.270 0.618 0.704 0.697 0.826 0.898...
         0.728 0.698 0.691 0.593 0.368 0.149 0.410 0.667 0.837...
         0.668 0.427 0.066 0.324 0.631 0.688 0.636 0.764 0.909...
         0.820 0.669 0.705 0.597 0.327 0.178 0.278 0.711 0.849];
B_exp = [0.028 0.141 0.231 0.292 0.018 0.020 0.041 0.252 0.098 0.052...
         0.087 0.133 0.139 0.059 0.152 0.273 0.176 0.068 0.027...
         0.049 0.194 0.153 0.047 0.070 0.096 0.196 0.077 0.018...
         0.041 0.146 0.061 0.013 0.068 0.229 0.198 0.222 0.002];
C_exp = [0.021 0.069 0.291 0.429 0.399 0.103 0.114 0.037 0.021 0.021...
         0.026 0.039 0.038 0.164 0.376 0.494 0.316 0.082 0.024...
         0.120 0.268 0.506 0.398 0.200 0.044 0.027 0.052 0.002...
         0.055 0.033 0.043 0.297 0.425 0.479 0.322 0.083 0.031];
     
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

ylabel({'p(1\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

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

A_exp = [0.914 0.668 0.208 0.069 0.022 0.070 0.230 0.545 0.670 0.918...
         0.771 0.452 0.389 0.121 0.094 0.048 0.203 0.703 0.862 0.727...
         0.385 0.023 0.054 0.184 0.376 0.594 0.784 0.906 0.812 0.513...
         0.397 0.155 0.057 0.040 0.242 0.623 0.894];     
B_exp = [0.006 0.036 0.113 0.031 0.080 0.283 0.486 0.351 0.121 0.031...
         0.056 0.318 0.440 0.377 0.118 0.051 0.128 0.038 0.044 0.133...
         0.445 0.685 0.844 0.690 0.639 0.240 0.069 0.060 0.070 0.262...
         0.544 0.731 0.834 0.812 0.389 0.214 0.025];  
C_exp = [0.035 0.202 0.413 0.779 0.759 0.410 0.065 0.006 0.058 0.042...
         0.040 0.012 0.070 0.456 0.784 0.757 0.326 0.113 0.010 0.029...
         0.085 0.061 0.103 0.056 0.035 0.040 0.060 0.012 0.044 0.010...
         0.027 0.046 0.034 0.111 0.098 0.107 0.058];     
     
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

ylabel({'p(1\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

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


A_exp = [0.856 0.687 0.330 0.029 0.076 0.161 0.315 0.519 0.739 0.873...
         0.679 0.508 0.334 0.094 0.041 0.040 0.353 0.672 0.878 0.688...
         0.346 0.061 0.061 0.177 0.358 0.473 0.780 0.901 0.788 0.503...
         0.321 0.087 0.006 0.066 0.319 0.734 0.872];
B_exp = [0.018 0.176 0.434 0.790 0.872 0.727 0.523 0.419 0.144 0.027...
         0.078 0.293 0.523 0.720 0.866 0.652 0.392 0.182 0.044 0.061...
         0.084 0.036 0.144 0.298 0.466 0.245 0.085 0.006 0.075 0.384...
         0.410 0.243 0.040 0.027 0.037 0.029 0.040];     
C_exp = [0.038 0.060 0.187 0.159 0.047 0.003 0.013 0.011 0.036 0.036...
         0.015 0.049 0.028 0.035 0.057 0.081 0.082 0.034 0.012 0.098...
         0.404 0.853 0.762 0.366 0.092 0.006 0.011 0.012 0.031 0.033...
         0.101 0.408 0.792 0.893 0.618 0.153 0.035];        
     
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

ylabel({'p(1\rightarrow n)'},'FontWeight','normal','FontSize',fontsize,...
    'FontName','Calibri');
set(gca,'Xlim',[0 2*pi],'Ylim',[0 1],'XTick',[0:pi/2:2*pi],'XTickLabel',{'0','0.5¦Ð','¦Ð','1.5¦Ð','2¦Ð'},'FontName','Calibri','FontSize',fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%