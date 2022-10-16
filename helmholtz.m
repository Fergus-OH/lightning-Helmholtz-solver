function [u, maxerr, tsolve, nkv, Z, Zplot, pol, A] = helmholtz(wavenum, P, varargin)
%HELMHOLTZ  Lightning Helmholtz solver.
%         U = HELMHOLTZ(wavenum,P,G) solves the Helmholtz equation with
%         Dirichlet boundary data on the simply-connected region Omega
%         bounded by P, which may be a polygon or circular polygon,
%         v1, University of Oxford, April 2022.
% 
%  Inputs:
%      P = vector of corners as complex numbers z = x+iy in counterclockwise
%              order to specify a polygon
%          or cell array of corners v and pairs [v r] to specify a circular
%              polygon: r = radius of curvature of arc from this v to the next
%          or one of the following specified strings
%             'sqr'[square], 'rec'[tangle], 'snow'[flake], pent[agaon],
%              'hex'[agon], 'L', 'circleL', or 'C'
%          or integer >= 3, the number of corners of a random polygon
%          or integer <= -3, -1 x no. of corners of a random circular polygon]
% 
%      g = function handle for Dirichlet boundary data that satisfies helm(g)
%          or cell array of function handles for sides P1-P2, P2-P3,...
%          (default @(z) exp(-1i*real(wavenum*exp(-1i*z0ang)*z))) for wavenumber>0
%                 @(z) @(z) besselh(0,-wavenum*abs(z-(z0_pt))) for wavenumber<0)
%
%  Further optional inputs:
%    'tol', tol = tolerance for maximal absolute error (default 1e-6)
%    'z0' to specify the point/angle of incidence for defulat incident source
%    'noplots' to suppress plotting
%    'noplot3d' to suppress 3D surface plotting
%    'steps' for step-by-step plots of errors on boundary and poles
%    'scat' to plot only the scattered field
%    'slow' to turn off adaptive mode for cleaner root-exp convergence curves
%    'fs' to set font size for plots
% 
%  Outputs for [U,MAXERR,TSOLVE,NKV,Z,ZPLOT,POL,A] = helmholtz(k,P,G)
%       u = function handle for solution u(z) of helm(u) = 0, u = g on boundary
%  maxerr = estimated upper bound on maximal error,
%             even near singularities at corners
%       tsolve = runtime for calculcating solution
%       nkv = vector containing number of poles at each corner in solution
%       pol = set of poles
%       Z = sample points on boundary
%       Zplot = set of points for plotting boundary
%       A = final rectangular matrix used for least-squares problem
% 
% Examples:
%
%   helmholtz(50,'sqr');                              % plane with square
%   helmholtz(-50,'sqr');                             % point with square
%   helmholtz(-20,'pent','tol',1e-10);                % pentagon
%   helmholtz(-20,'circleL','z0',2+3i);               % circular L-shape
%   helmholtz(20,'bullet','z0',1);                    % bullet
%   helmholtz(-30,'snow','steps')                     % snowflake
%   helmholtz(50,[1/2*exp(2i*pi*([1:3])/3)],'z0',1i)  % triangle
% 
% two point sources:
%   wavenum = -30; z0_pt = .5+1i;
%   g = @(z) besselh(0,-wavenum*abs(z-(z0_pt))) + ...
%       besselh(0,-wavenum*abs(z-(-z0_pt')));
%   helmholtz(wavenum,'sqr',g,'noplot3d');

%% Set up the problem
[g, P, w, ww, pt, dw, tol, steps, scat, ...      % parse inputs
   plots, plot3d, fs, slow] = ...
   parseinputs(wavenum, P, varargin{:});
wavenum = abs(wavenum);                         
Zplot = ww;
nw = length(w);                                  % number of corners
x = real(ww); y = imag(ww);
wr = sort(x); wr = wr([1 end]);
wi = sort(y); wi = wi([1 end]);
scl = max([diff(wr),diff(wi)]);                  % characteristic length scale
q = .5; if slow == 1, q = 0; end                 % sets which corners get more poles
inpolygonc = @(z,w) inpolygon(real(z), ...       % complex variant of "inpolygon"
            imag(z),real(w),imag(w));  
                                                   
A = x(1:end).*y([2:end,1])-x([2:end,1]).*y(1:end); 
As = sum(A)/2;
x_bar = (sum((x([2:end,1])+x(1:end)).*A)*1/6)/As;
y_bar = (sum((y([2:end,1])+y(1:end)).*A)*1/6)/As;
wc = x_bar+1i*y_bar;                              % centroid of polygon
zs = wc;
while ~inpolygonc(zs,ww)
   zs = input('centroid not in polygon, please manually enter runge point: ');
end
for k = 1:nw
   forward = pt{k}(.01*dw(k)) - w(k);            % small step toward next corner
   j = mod(k-2,nw)+1;
   backward = pt{j}(.99*dw(j)) - w(k);           % small step toward last corner
   tmp = -1i*backward*sqrt(-forward/backward);
   inward(k) = tmp/abs(tmp);                    % inward direction from corner
end
warn = warning('off','MATLAB:rankDeficientMatrix');  % matrices are ill-conditioned

%% Set up for plots
if plots
   LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
   PO = 'position'; FW = 'fontweight'; NO = 'normal';
   npts = max(120,min(wavenum*8,250));
   sx = linspace(wr(1)-scl,wr(2)+scl,npts); sy = linspace(wi(1)-scl,wi(2)+scl,npts);
   [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
   ax = [wr(1:2); wi(1:2)] + scl*[-1 1 -1 1]';
   axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';
   axnarrow = [wr(1:2); wi(1:2)] + .5*scl*[-1 1 -1 1]';
end
if steps
   figure, subplot(1,2,1), shg, plot(w([1:end 1]),'k',LW,1)
   grid on, axis equal, axis(axnarrow)
end

%% Main loop: increase number of poles until convergence ==============================
Nvec = []; errvec = []; tic
errk = ones(nw,1);                               % max error near each corner
nkv = ones(nw,1);                                % no. of poles at each corner
maxstepno = 30; degfmax = 35^2; err0 = Inf;

for stepno = 1:maxstepno
   % Fix poles and sample pts on bndry.  Side k means side from corner k to k+1.
   Z = [];           % col vector of sample points on boundary
   G = [];           % col vector of boundary values at these points
   pol = [];         % row vector of poles of the rational approximation
   J = [];           % row vector of indices of which corner each pole belongs to
   d = [];           % row vector of distances from poles to their corners
   tt = cell(nw,1);  % cell array of distances of sample points along each side
   for k = 1:nw
      nk = nkv(k);                                      % no. of poles at this corner
%       normal clustering
%             sk = sqrt(1:nk) - sqrt(nk);
%             dk1 = exp(4*sk); dk = scl*dk1;              % stronger clustering near corner
%             dk = dk(dk>1e-15*scl);                      % remove poles too close to corner

      % stronger clustering (greased lightning)
      sk = sqrt(1:nk) - sqrt(nk);
      dk1 = exp(4*sk);
      dk = .5*scl*dk1/max(dk1);                         % stronger clustering near corner
      dk = dk(dk>1e-15*scl);                            % remove poles too close to corner
      polk = w(k) + inward(k)*dk;                       % poles near this corner
      ii = find(~inpolygonc(polk(dk>1e-12*scl),ww),1);  % work around inaccuracy
      if length(ii)>0                                   % don't allow poles in Omega
         dk = dk(1:ii-2); polk = polk(1:ii-2);
      end
      pol = [pol polk]; d = [d dk];
      dvec = [(1/3)*dk (2/3)*dk dk];                    % finer pts for bndry sampling
      nsidepts = sqrt(dw(k))*max(30,10*sqrt(nk));       % numbr additional pts along side
      if length(P{k}) == 2                              % add more bndry pts if side curved
         factor = 1+dw(k)/abs(P{k}(2));
         nsidepts = 1.5*nsidepts*factor;
      end
      tt{k} = [tt{k} dvec(dvec<dw(k)) ...               % add clustered pts near corner
         linspace(0,dw(k),nsidepts)];
      j = mod(k-2,nw)+1;                                % index of last corner
      tt{j} = [tt{j} dw(j)-dvec(dvec<dw(j)) ...         % likewise in other direction
         dw(j)-linspace(0,dw(j),nsidepts)];             
   end
   for k = 1:nw
      tt{k} = sort(unique(tt{k}(:)));
      tk = tt{k}; pk = pt{k};                           % abbrevations
      Z = [Z; pk(tk)];                                  % sample pts on side k
      G = [G; g{k}(pk(tk))];                            % boundary data at these pts
   end

   % Solve the Helmholtz problem
   n = 4*stepno;                                        % (degree of polynomial term) i.e number of runge terms
   Np = length(pol);                                    % number of newman terms

   M = size(Z,1);                                       % number of boundary samples
   Gn = G;
   A = cauchy(Z,zs,pol,n,wavenum);
   J = [zeros(1,2*n+1) J J];                            % corner for each col
   N = size(A,2);                                       % no. of cols = 2n+1+2Np CHECK WHY THIS MAY NOT BE WORKING
   Kj = zeros(M,1);
   for j = 1:M
      dd = abs(Z(j)-w);
      Kj(j) = find(dd==min(dd),1);                      % nearest corner to Zj
   end
   D = 1./vecnorm(A); nD = length(D); D = spdiags(D',0,nD,nD);
   c = D*((A*D)\Gn);
   ftemp = @(zz) cauchy(zz,zs,pol,n,wavenum)*c;  
   u = @(z) reshape(ftemp(z(:)),size(z));               % to u and f both allowed
   for k = 1:nw
      Kk = find(Kj==k);
      errk(k) = norm(A(Kk,:)*c-Gn(Kk),inf);             % error near corner k
   end
   err = norm(A*c-Gn,inf);                              % global error
   polmax = 120;
   for k = 1:nw
      if (errk(k) > q*err) & (nkv(k) < polmax)
         nkv(k) = nkv(k)+ceil(1+sqrt(nkv(k)));          % increase no. poles
      else
         nkv(k) = max(nkv(k),ceil(stepno/2));
         %nkv(k) = min(polmax,nkv(k)+1);
      end
   end
   if steps                                             % plot error on bndry
      subplot(1,2,1), plot(ww,'k',LW,1), grid on
      title(sprintf('step %d',stepno))
      axis equal, axis(axnarrow), hold on
      plot(pol,'.r',MS,7), hold off
      subplot(1,2,2), semilogy([-pi pi],err*[1 1],'--b',LW,1), axis square
      hold on, axis([-pi pi 1e-16 100]), grid on
      semilogy(angle(Z-wc),abs(u(Z)-G),'.k',MS,4)
      semilogy(angle(pol-wc),d,'.r',MS,7), hold off
      set(gca,'ytick',10.^(-16:4:0))
      set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
      set(gca,FS,fs-1), xlabel('angle on boundary wrt wc',FS,fs)
      title('bndry err (black) & poles (red)',FS,fs)
      disp('Press <enter> for next plot'); pause
   end
   errvec = [errvec err]; Nvec = [Nvec; N];
   if err < tol, break, end                             % convergence success
   if err < err0                                        % save the best so far
      u0 = u; Z0 = Z; G0 = G; A0 = A; M0 = M;
      N0 = N; err0 = err; pol0 = pol;
   end
   if (N > degfmax) || (stepno == maxstepno) || (Np >= polmax*nw)  % failure
      u = u0; Z = Z0; G = G0; A = A0; M = M0;
      N = N0; err = err0; pol = pol0;
      warning('HELMHOLTZ failure.  Loosen tolerance or add corners?')
      break
   end
end
warning(warn.state,'MATLAB:rankDeficientMatrix')       % back to original state
tsolve = toc;  % =========== end of main loop =========================================

%% Finer mesh for a posteriori error check
Z2 = []; G2 = [];
for k = 1:nw
   newtt = mean([tt{k}(1:end-1) tt{k}(2:end)],2);
   newpts = pt{k}(newtt);
   Z2 = [Z2; newpts];
   G2 = [G2; g{k}(newpts)];
end
err2 = norm((G2-u(Z2)),inf);
maxerr = max(err,err2);                                % estimated max error

%% Convergence curve plot
if plots
   fig = figure; fig.Position = [500 100 700 400]; shg
   axes(PO,[.07 .7 .3 .2])
   semilogy(sqrt(Nvec),errvec,'.-k',LW,0.7,MS,10), grid on, hold on
   semilogy(sqrt(N),maxerr,'or',MS,7,LW,1), hold off
   errmin = .01*tol; axis([0 1.1*max(sqrt(Nvec)) 1e-14 100])

   set(gca,FS,fs-1), title('convergence',FS,fs,FW,NO)
   xlabel('sqrt(DoF)',FS,fs), ylabel('error',FS,fs)
   set(gca,'ytick',10.^(-16:4:0))
   ax2 = axis; x1 = ax2(1) + .05*diff(ax2(1:2));
   s = sprintf('solve time = %6.3f secs',tsolve);

   if ~steps, text(x1,4e-11,s,FS,fs), end
   z = randn(1000,1)+1i*randn(1000,1); z = z/10;
   tic, u(z); teval = 1e3*toc;
   s = sprintf('eval time = %4.1f microsecs per pt',teval);
   text(x1,4e-13,s,FS,fs)
end

%% Error plot along boundary
if plots
   axes(PO,[.07 .4 .3 .2])
   semilogy([-pi pi],maxerr*[1 1],'--b',LW,1), hold on
   semilogy(angle(Z2-wc),abs(u(Z2)-G2),'.r',MS,4)
   axis([-pi pi .0001*errmin 1]), grid on
   semilogy(angle(Z-wc),abs(u(Z)-G),'.k',MS,4), hold off
   set(gca,'ytick',10.^(-16:4:0))
   set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
   set(gca,FS,fs-1), xlabel('angle on boundary wrt wc',FS,fs)
   ylabel('error',FS,fs)
   title('error on boundary',FS,fs,FW,NO)
   %    axis square
end

%% Displaying number of poles per corner
if plots
   axes(PO,[.07 .1 .3 .2])
   bar(angle(w-wc),nkv)
   set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
   set(gca,FS,fs-1)
   xlabel('angle on boundary wrt wc',FS,fs)
   ylabel('#poles',FS,fs)
   xlim([-pi,pi])
   ylim([0 round(max(nkv),-1)+10])
   title('pole distribution',FS,fs,FW,NO)
end

%% Plot solution of total field
if plots
   tic,
   uu = u(zz(:));
   ind = inpolygonc(zz,ww);
   mxcol = .5*max(real(uu(~ind))); mncol = min(real(uu(~ind)));
   uu = reshape(uu,size(zz));
   if scat, ff = g{1}(zz) - uu;
   else ff = uu; end
   axes(PO,[.40 .1 .6 .8])
   imagesc(sx,sy,real(exp(-1i*0)*ff),[mncol mxcol]), hold on
   fill(real(ww),imag(ww),[1 1 1]), hold on, colorbar, set(gca,FS,8)
   set(gca,'YDir','normal')
   mxabscol = max(abs([mxcol,mncol]));
   caxis([-mxabscol,mxabscol])
   axis equal, axis(ax)

   plot(pol,'.r',MS,8), plot(real(zs),imag(zs),'.k',MS,8), plot(Z,'.k',MS,2), hold off
   ss = ['wavenum = ' int2str(wavenum) ', dim(A) = ' int2str(M) ' x ' int2str(N) ', ', ...
      ' #poles = ' int2str(length(pol))];
   title(ss,FS,fs,FW,NO)
end

%% 3D Plotting
if plots && plot3d
   sz = real(exp(1i*0)*ff); sz(ind) = nan;
   figure, axis equal, surfc(sx,sy,sz)
end

end
%% Functions
%% Treat the inputs
function [g, P, w, ww, pt, dw, tol, steps, scat, ...
   plots, plot3d, fs, slow] = parseinputs(wavenum, P, varargin)
%% Defaults
tol = 1e-6; steps = 0; scat = 1; plots = 1;
plot3d = 1; fs = 9; slow = 0;

z0_pt = .5+1i; z0ang = 5*pi/6;

%% First treat the domain, defined by P
randomcirc = 0;
% this is if P is not a cell i.e P!={[..]}
if ~iscell(P)
   if isnumeric(P)
      if length(P) > 1, w = P;                    % vertices have been specified
      else
         if P < 0
            randomcirc = 1; P = -P;               % random circular arcs
         end
         w = exp(2i*pi*(1:P)/P).*(.1+rand(1,P));  % random vertices
      end
   else
      cellscale = @(P,s) cellfun(@times, P, num2cell(s*ones(1,length(P))), 'UniformOutput',false);
      if strcmp(P,'sqr'), w = .5*[-1-1i 1-1i 1+1i -1+1i];
      elseif strcmp(P,'rec'), w = .5*[-2-1i 2-1i 2+1i -2+1i];
      elseif strcmp(P,'snow'), P = exp(2i*pi*(1:12)/12);
         w = P.*(1+.2*(-1).^(1:12)); w = w/2.8;
      elseif strcmp(P,'pent'), w = .7*exp(pi*2i*(1:5)/5);
      elseif strcmp(P,'hex'), w = .7*exp(pi*2i*(1:6)/6);
      elseif strcmp(P,'kite'), w = .25*[0 2+4i 5i -2+4i];
      elseif strcmp(P,'L'), w = .5*([2 2+1i 1+1i 1+2i 2i 0]-.5*(1+1i));
      elseif strcmp(P,'circleL'), w = {2 [2+1i -1] 1+2i 2i 0};
      elseif strcmp(P,'C'), P = {-2-1i 2-1i 2+1i [1+1i -1.1] -1+1i -2+1i}; P = cellscale(P,.25);
      elseif strcmp(P, 'bullet'), P = {.5*[1-.5i .6] .5*(1+.5i) .5*(-1+.5i) .5*(-1-.5i)};
      end
   end
   if ~iscell(P), P = num2cell(w); end            % convert to cell array
   if randomcirc
      for k = 1:length(P)
         r = .6/rand;
         P{k} = [P{k} r*(-1)^double(randn>0)];
      end
   end
end

nw = length(P);
for k = 1:nw
   w(k) = P{k}(1);
end
w = w(:);
ww = [];                                            % bndry pts for plotting
for k = 1:nw
   kn = mod(k,nw)+1;                                % index of next corner
   ww = [ww; w(k)];
   if isnumeric(P{k})
      if length(P{k}) == 1                          %     straight arc
         dw(k) = abs(w(kn)-w(k));                   % distance to next corner
         pt{k} = @(t) w(k) + t*(w(kn)-w(k))/dw(k);  % parametrization of arc
      else                                          %     circular arc
         r = P{k}(2);                               % radius of arc
         a = w(k); b = w(kn); ab = abs(b-a);        % endpoints of arc
         theta = real(asin(ab/(2*r)));              % half-angle of arc
         c = a + r*exp(1i*(pi/2-theta))*(b-a)/ab;   % center of arc
         dw(k) = 2*theta*r;                         % arc length of arc
         pt{k} = @(t) c - ...
            r*exp(1i*(pi/2+t/r-theta))*(b-a)/ab;   % parametrization of arc
         ww = [ww; pt{k}(linspace(0,dw(k),50)')];
      end
   else
      error('HELMHOLTZ:parseinputs','general boundary arcs not yet implemented')
   end
end
ww = [ww; w(1)];
Zplot = ww;

%% Next treat the boundary conditions
for k = 1:nw
   % defaults
   if wavenum>0, g{k} = @(z) exp(-1i*real(wavenum*exp(-1i*z0ang)*z));
   else g{k} = @(z) besselh(0,-wavenum*abs(z-(z0_pt))); end
end

j = 1;
while j < nargin-1
   j = j+1;
   v = varargin{j-1};

   if ~ischar(v)                 % This block specifies Dirichlet bndry data g.
      if isa(v,'cell')           % if cell array, nothing to change
         g = v;
      elseif isa(v,'double')     % if vector, convert to cell array of fun. handles
         for k = 1:nw
            g{k} = @(z) v(k) + 0*z;
         end
      elseif isa(v,'function_handle')  % if fun. handle, convert to cell array
         for k = 1:nw
            g{k} = @(z) v(z);
         end
      else
         error('HELMHOLTZ:parseinputs','boundary data g not in correct form')
      end
      
   elseif strcmp(v,'z0'), j = j+1;  % user-specified z0 parameter for default bndry funs
      z0param = varargin{j-1};
      for k = 1:nw
         if wavenum>0
            if ~isreal(z0param), z0ang = angle(z0param); else z0ang = z0param; end
            g{k} = @(z) exp(-1i*real(wavenum*exp(-1i*z0ang)*z));
         else z0_pt = z0param; g{k} = @(z) besselh(0,-wavenum*abs(z-(z0_pt))); end
      end

   elseif strcmp(v,'tol'), j = j+1; tol = varargin{j-1};
   elseif strcmp(v,'steps'), steps = 1; plots = 1;
   elseif strcmp(v,'noscat'), scat = 0;
   elseif strcmp(v,'noplots'), plots = 0;
   elseif strcmp(v,'noplot3d'), plot3d = 0;
   elseif strcmp(v,'slow'), slow = 1;
   elseif strcmp(v,'fs'), j = j+1; fs = varargin{j-1};
   else error('HELMHOLTZ:parseinputs','Unrecognized string input')
   end
end
wavenum = abs(wavenum);

end   % end of parseinputs

%% Cauchy matrix
function C = cauchy(Z,zs,pol,n,wavenum)            % Make Cauchy matrix
Np = length(pol);
M = length(Z);
C = zeros(M,2*(Np+n)+1);
if Np>0                                            % Newman parts
   D = Z-pol; absD = abs(D);
   vals = [real(D./absD) imag(D./absD)];
   bessvals = besselh(1,wavenum*absD);
   C(:,1:2*Np) = repmat(bessvals,[1,2]).*vals;
end

s = Z-zs; abss = abs(s);
for j = 1:n                                        % Runge parts
   bessval = besselh(j,wavenum*abss);
   C(:,2*Np+j) = bessval.*real((s./abss).^j);
   C(:,2*Np+n+j) = bessval.*imag((s./abss).^j);
end
C(:,2*(Np+n)+1) = besselh(0,wavenum*abss);

end