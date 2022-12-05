function [Q,fcnt] = myquad(funfcn,a,b,tol,trace,varargin)
%QUAD   Numerically evaluate integral, adaptive Simpson quadrature.
%   Q = QUAD(FUN,A,B) tries to approximate the integral of function
%   FUN from A to B to within an error of 1.e-6 using recursive
%   adaptive Simpson quadrature.  The function Y = FUN(X) should
%   accept a vector argument X and return a vector result Y, the
%   integrand evaluated at each element of X.  
%
%   Q = QUAD(FUN,A,B,TOL) uses an absolute error tolerance of TOL 
%   instead of the default, which is 1.e-6.  Larger values of TOL
%   result in fewer function evaluations and faster computation,
%   but less accurate results.  The QUAD function in MATLAB 5.3 used
%   a less reliable algorithm and a default tolerance of 1.e-3.
%
%   [Q,FCNT] = QUAD(...) returns the number of function evaluations.
%
%   QUAD(FUN,A,B,TOL,TRACE) with non-zero TRACE shows the values
%   of [fcnt a b-a Q] during the recursion.
%
%   QUAD(FUN,A,B,TOL,TRACE,P1,P2,...) provides for additional 
%   arguments P1, P2, ... to be passed directly to function FUN,
%   FUN(X,P1,P2,...).  Pass empty matrices for TOL or TRACE to
%   use the default values.
%
%   Use array operators .*, ./ and .^ in the definition of FUN
%   so that it can be evaluated with a vector argument.
%
%   Function QUADL may be more efficient with high accuracies
%   and smooth integrands.
%
%   Example:
%       FUN can be specified three different ways.
%
%       A string expression involving a single variable:
%          Q = quad('1./(x.^3-2*x-5)',0,2);
%
%       An inline object:
%          F = inline('1./(x.^3-2*x-5)');
%          Q = quad(F,0,2);
%
%       A function handle:
%          Q = quad(@myfun,0,2);
%          where myfun.m is an M-file:
%             function y = myfun(x)
%             y = 1./(x.^3-2*x-5);
%
%   See also QUADL, DBLQUAD, INLINE, @.

%   Based on "adaptsim" by Walter Gander.  
%   Ref: W. Gander and W. Gautschi, "Adaptive Quadrature Revisited", 1998.
%   http://www.inf.ethz.ch/personal/gander
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 1.1.1.1 $  $Date: 2002/12/05 20:52:54 $

f = fcnchk(funfcn);
if nargin < 4 | isempty(tol), tol = 1.e-6; end;
if nargin < 5 | isempty(trace), trace = 0; end;

% Initialize with three unequal subintervals.
h = 0.13579*(b-a);
x = [a a+h a+2*h (a+b)/2 b-2*h b-h b];
y1 = feval(f, x(1), varargin{:});
y2 = feval(f, x(2), varargin{:});
y3 = feval(f, x(3), varargin{:});
y4 = feval(f, x(4), varargin{:});
y5 = feval(f, x(5), varargin{:});
y6 = feval(f, x(6), varargin{:});
y7 = feval(f, x(7), varargin{:});
fcnt = 7;

% Fudge endpoints to avoid infinities.
if ~min(min(min(isfinite(y1))))
   y1 = feval(f,a+eps*(b-a),varargin{:});
   fcnt = fcnt+1;
end
if ~min(min(min(isfinite(y7))))
   y7 = feval(f,b-eps*(b-a),varargin{:});
   fcnt = fcnt+1;
end
Q=zeros(3,size(y1,1),size(y1,2),size(y1,3));
% Call the recursive core integrator.
hmin = eps/1024*abs(b-a);
[Q(1,:,:,:),fcnt,warn(1)] = ...
   quadstep(f,x(1),x(3),squeeze(y1),squeeze(y2),squeeze(y3),tol,trace,fcnt,hmin,varargin{:});
[Q(2,:,:,:),fcnt,warn(2)] = ...
   quadstep(f,x(3),x(5),squeeze(y3),squeeze(y4),squeeze(y5),tol,trace,fcnt,hmin,varargin{:});
[Q(3,:,:,:),fcnt,warn(3)] = ...
   quadstep(f,x(5),x(7),squeeze(y5),squeeze(y6),squeeze(y7),tol,trace,fcnt,hmin,varargin{:});
Q = squeeze(sum(Q,1));
warn = max(warn);

switch warn
   case 1
      warning('Minimum step size reached; singularity possible.')
   case 2
      warning('Maximum function count exceeded; singularity likely.')
   case 3
      warning('Infinite or Not-a-Number function value encountered.')
   otherwise
      % No warning.
end

% ------------------------------------------------------------------------

function [Q,fcnt,warn] = quadstep (f,a,b,fa,fc,fb,tol,trace,fcnt,hmin,varargin)
%QUADSTEP  Recursive core routine for function QUAD.

maxfcnt = 10000;

% Evaluate integrand twice in interior of subinterval [a,b].
h = b - a;
c = (a + b)/2;
if abs(h) < hmin | c == a | c == b
   % Minimum step size reached; singularity possible.
   Q = h*fc;
   warn = 1;
   return
end
x = [(a + c)/2, (c + b)/2];
y1 = feval(f, x(1), varargin{:});
y2 = feval(f, x(2), varargin{:});

fcnt = fcnt + 2;
if fcnt > maxfcnt
   % Maximum function count exceeded; singularity likely.
   Q = h*fc;
   warn = 2;
   return
end
fd = y1;
fe = y2;

% Three point Simpson's rule.
Q1 = (h/6)*(fa + 4*fc + fb);

% Five point double Simpson's rule.
Q2 = (h/12)*(fa + 4*fd + 2*fc + 4*fe + fb);

% One step of Romberg extrapolation.
Q = Q2 + (Q2 - Q1)/15;

if ~min(min(min(isfinite(Q))))
   % Infinite or Not-a-Number function value encountered.
   warn = 3;
   return
end
if trace
   disp(sprintf('%8.0f %16.10f %18.8e %16.10f',fcnt,a,h,Q))
end

% Check accuracy of integral over this subinterval.

if max(max((abs(Q2 - Q)))) <= tol
   warn = 0;
   return

% Subdivide into two subintervals.
else
   [Qac,fcnt,warnac] = quadstep(f,a,c,fa,fd,fc,tol,trace,fcnt,hmin,varargin{:});
   [Qcb,fcnt,warncb] = quadstep(f,c,b,fc,fe,fb,tol,trace,fcnt,hmin,varargin{:});
   Q = Qac + Qcb;
   warn = max(warnac,warncb);
end
