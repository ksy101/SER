function [x,u,res_p,res_d] = reg_TV2(A,y,T,varargin)

% x and y directions gradient

[LM,n] = size(A);
% data set size
[B,p] = size(y);
if (LM ~= B)
    error('mixing matrix M and data set y are inconsistent');
end


%%
% maximum number of AL iteration
AL_iters = 300;
% regularizatio parameter
lambda = 4e-2;
% Positivity constraint
positivity = 'yes';
% tolerance for the primal and dual residues
tol = 1e-6;
% initialization
x0 = 0;
rho = 1;
%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'LAMBDA'
                lambda = varargin{i+1};
                if (sum(sum(lambda < 0)) >  0 )
                       error('lambda must be positive');
                end
            case 'RHO'
                rho = varargin{i+1};
                if (sum(sum(rho < 0)) >  0 )
                       error('rho must be positive');
                end                
            case 'POSITIVITY'
                positivity = varargin{i+1};
            case 'TOL'
                tol = varargin{i+1};
            case 'X0'
                x0 = varargin{i+1};
                if (size(x0,1) ~= n) | (size(x0,1) ~= p)
                    error('initial X is  inconsistent with M or Y');
                end
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

%---------------------------------------------
%  If lambda is scalar convert it into vector
%---------------------------------------------
Nlambda = size(lambda);
if Nlambda == 1
    % same lambda for all pixels
    lambda = lambda*ones(n,p);
elseif Nlambda ~= p
        error('Lambda size is inconsistent with the size of the data set');
else
  %each pixel has its own lambda
   lambda = repmat(lambda(:)',n,1);
end

% compute mean norm
  norm_y = sqrt(mean(mean(y.^2)));
% rescale M and Y and lambda
if norm_y ~=0
    A = A/norm_y;
    y = y/norm_y;
    lambda = lambda/norm_y^2;
else
    A=A;
    y=y;
    lambda=lambda;
end

  


%%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------
mu_AL = 0.01;
muu = 10*mean(lambda(:)) + mu_AL;

% rho = 1;
F = A'*A+rho*eye(n);

IF = pinv(F);

e = ones(p,1);
D = spdiags([e -e], -1:0, p, p);
np = sqrt(p);
for i =1:np
    for j=1:np
        D(np*i,np*j)=0;
        D(min(np*i+1,p),np*j)=0;
    end
end
    Td=T*D;

Ix = pinv(eye(p)+D*D'+(Td)*(Td)');

%%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
if x0 == 0
    x= IF*A'*y;
end

u = x;
% dx = gradient(x);
dx = x*D;
% dy = gradient(x*T);
dy = x*Td;
% scaled Lagrange Multipliers
zu  = 0*u;
zx  = 0*dx;
zy  = 0*dy;


%%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt(p*n)*tol;
tol2 = sqrt(p*n)*tol;
i=1;
res_p = inf;
res_d = inf;
maskx = ones(size(x));
mu_changed = 0;
while (i <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2)) 
    if mod(i,10) == 1
        dx0 = dx;
        dy0 = dy;
    end
    
    u = IF*(A'*y+rho*(x+zu));
    gx = x*D;
    dx =  soft(gx+zx,lambda/rho);
    gy = x*Td;
    dy =  soft(gy+zy,lambda/rho);
    
    x = ((u-zu)+(dx-zx)*D'+(dy-zy)*(Td)')*Ix;

    if strcmp(positivity,'yes')
       maskx = (x >= 0);
       x = x.*maskx; 
%        xw=x(1,:);
%        xw(xw>1)=1;
%        x(1,:) = xw;
%        xp=x(2,:);
%        xp(xp>1.18)=1.18;
%        x(2,:) = xp;   
    end
    
    zu = zu + x - u; 
    zx = zx +dx - gx;
    zy = zy +dy - gy;

    res_p = norm(x-u,'fro');
    res_d = muu*norm(dx-dx0,'fro');
    
%     res_gt = norm(gt-x);
%     disp([strcat('In iteration ',' ', num2str(i)),';','error is',':',num2str(res_gt)]);
    
    i=i+1;

       
end

    
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
