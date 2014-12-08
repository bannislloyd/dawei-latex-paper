function[sigma_x,sigma_y,sigma_z,I,Ix,Iy,Iz,ST0,ST1,KIx,KIy,KIz] = MultiPauli(n)

sigma_x=[0,1;1,0];                  %sigma算符 
sigma_y=[0,-i;i,0];
sigma_z=[1,0;0,-1];
ST0=[1,0;0,0];
ST1=[0,0;0,1];
Ix=0.5*sigma_x;                      %单核的核自旋算符
Iy=0.5*sigma_y;
Iz=0.5*sigma_z;
I=[1,0;0,1];

%-----------华丽的分割线------------

KIx = cell(n,1);

for ii = 1:n
    KIx{ii} = 1;
end


for ii = 1:n
    
    for jj = 1:n
       if jj==ii
          KIx{ii} = kron(KIx{ii},Ix);
       else 
          KIx{ii} = kron(KIx{ii},I);
       end
    end
end


KIy = cell(n,1);

for ii = 1:n
    KIy{ii} = 1;
end


for ii = 1:n
    
    for jj = 1:n
       if jj==ii
          KIy{ii} = kron(KIy{ii},Iy);
       else 
          KIy{ii} = kron(KIy{ii},I);
       end
    end
end        

KIz = cell(n,1);

for ii = 1:n
    KIz{ii} = 1;
end


for ii = 1:n
    
    for jj = 1:n
       if jj==ii
          KIz{ii} = kron(KIz{ii},Iz);
       else 
          KIz{ii} = kron(KIz{ii},I);
       end
    end
end