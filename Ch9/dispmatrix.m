function dm = dispmatrix( sites, type, varargin )
% DISPMATRIX This function generates various types of dispersal matrices.

%
% Currently supported: gaussian, exponential and identity.
%
% Usage:
%
% dm = dispmatrix (numsites, type, parameter(s), shift, boundaries)
% 
% where type can be 'gaussian', 'exponential', 'pool', 'twobay', 'uniform'
% or 'identity'.  The other inputs are not necessary in all cases.
% Parameter represents the relevant parameter in the type used - sigma for
% gaussian or exponent for exponential.  In twobay, several parameters
% are needed.
%
% The pool case corresponds to larvae evenly distributed between all
% cases.  The identity case corresponds to no dispersal.
%
% NOTE: boundaries is a string which tells what to do at the ends.  It
% can be one of the following: 'nothing' - just let things fall off end;
% 'conservative' - which normalizes the 'nothing' result to all larvae
% land in system; and 'circular' - the last site goes to the first site.
% The default is nothing.
%
% Shift gives the shift from center and is assumed zero if not given.
%

% Look through type for various possibilities.
switch type
 case 'gaussian'
  % Gaussian dispersal matrix with a possible shift.
  if nargin < 3, error('DMK: Too few arguments.'); end
  sigma = varargin{1};
  shift = 0;
  boundaries = 'nothing';
  
  if nargin > 3, shift = varargin{2}; end
  if nargin > 4, boundaries = varargin{3}; end

  % if sigma & shift are not scalars (i.e., non-uniform case), need to make adjustment
  if length(sigma) > 1
      sigma = repmat(sigma(:)',[sites,1]);
  end
  
  if length(shift) > 1
      shift = repmat(shift(:)',[sites,1]);
  end
  
  
  dm = gaussdisp ( sites, sigma, shift, boundaries );
  
 case { 'exponential', 'laplace', 'laplacian' }
  % This is the standard negative exponential case.
  % Can't remember what it is really called.
  if nargin < 3, error('DMK: Too few arguments.'); end
  a_exp = varargin{1};
  shift = 0;
  boundaries = 'nothing';
  
  if nargin > 3, shift = varargin{2}; end
  if nargin > 4, boundaries = varargin{3}; end

  dm = exponentialdisp ( sites, a_exp, shift, boundaries );
  
 case { 'pool' , 'LPER' }
  % Case where larvae evenly distributed among all sites.
  dm = ones( sites ) / sites;
  
 case { 'none' , 'identity' }
  % Just the identity matrix.  No dispersal case.
  dm = eye ( sites );

 case 'twobay'
  % Two embayment dispersal matrix.
  if nargin < 4, error('DMK: Too few arguments.'); end
  c1 = varargin{1};
  c2 = varargin{2};
  N = round( sites / 2 );
  leak = 0;
  
  if nargin > 4, N = varargin{3}; end
  if nargin > 5, leak = varargin{4}; end
  if N > sites, error('DMK: Second embayment beyond site limit.'); end
  
  dm = twobaydisp( sites, N, c1, c2, leak );
  
 case { 'even', 'uniform' }
  % This is dispersal uniformly within a distance.
  if nargin < 3, error('DMK: Too few arguments.'); end
  shift = 0;
  boundaries = 'nothing';
  
  dist = varargin{1};
  if nargin > 3, shift = varargin{2}; end
  if nargin > 4, boundaries = varargin{3}; end
  
  dm = uniformdisp( sites, dist, shift, boundaries );
 
 otherwise
  error('DMK: Unknown type of dispersal matrix.');
end


function dm = exponentialdisp ( sites, a_exp, shift, boundaries )
% EXPONENTIALDISP Generate a negative exponential dispersal matrix.
%
% Usage:
%
% dm = exponentialdisp (numsites, a, shift, boundaries)
% 

if nargin < 4
  boundaries = 'nothing';
end

if nargin < 3
  shift = 0;
end

x = repmat(1:sites,[sites,1]);

x = x' - x - shift;
  
dm = zeros( sites );

% Different types of boundary conditions.
switch boundaries
 case 'conservative' % Normalize larvae.  
  % The 1 / a_exp fixes fact that it should be 1/2a * exp(-x/a)
  dm = expint ( 1 / a_exp, x-0.5, x+0.5 );

  dm = dm ./ repmat(sum(dm),[sites,1]);
 case 'circular' % Loop boundary conditions.
  % Go out to 10 loops in either direction. 
  for i = -10:10
    dm = dm + expint( 1 / a_exp, x+i*sites-0.5, x+i*sites+0.5 );
  end
  
  % circular implies that it should be conservative.
  dm = dm ./ repmat(sum(dm),[sites,1]);
 case 'nothing'
  % The 1 / a_exp fixes fact that it should be 1/2a * exp(-x/a)
  dm = expint ( 1 / a_exp, x-0.5, x+0.5 );
 otherwise
  % The 1 / a_exp fixes fact that it should be 1/2a * exp(-x/a)
  dm = expint ( 1 / a_exp, x-0.5, x+0.5 );

  % Don't know what to do.
  warning( 'DMK - unknown boundary condition.' );
end

% This function integrates a_exp/2 * exp(-a_exp * abs(x)) from start to
% end.
%
% Note that the end must be greater than the start
function d = expint ( a_exp, s, e )
d = zeros(size(s));

% Deal with case of zero dispersal distance.
if isinf( a_exp )
  if a_exp > 0
    a_exp = realmax;
  else
    a_exp = realmin;
  end
end

posS = s>=0;
posE = e>=0;

pos = posS & posE;
neg = (~posS) & (~posE);
mixed = (~pos) & (~neg);

% Have to go by case due to absolute value in function.
d(pos)   = 0.5 * ( exp(-a_exp * s(pos)) - ...
			   exp(-a_exp * e(pos)) );

d(neg) = 0.5 * ( exp( a_exp * e(neg)) - ...
			 exp( a_exp * s(neg)) );

% Mixed a bit harder involving two integrals.
d(mixed) = 1 - 0.5 * ( exp(a_exp * s(mixed)) + exp(-a_exp * e(mixed)) );

function p = laplacian( x, mu, sigma )
%function p = laplacian( x, mu, sigma )

p = exp( -abs(x-mu)./sigma ) ./ 2 ./ sigma;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dm = gaussdisp ( sites, sigma, shift, boundaries )

% generate gaussian dispersal matrix


if nargin < 4
  boundaries = 'nothing';
end

if nargin < 3
  shift = 0;
end

sigma = max(sigma,realmin);

x = repmat(1:sites,[sites,1]);

x = x' - x;
  
dm = zeros( sites );

switch boundaries
    case 'nothing'
        
        dm = gaussdm(sigma, shift, x-0.5, x+0.5);
        
    case 'conservative'
        
        dm = gaussdm(sigma, shift, x-0.5, x+0.5);
        dm = dm./repmat(sum(dm),[length(dm),1]);
        
    case 'circular'
        
        for i = -10:1:10 % loop over 10 shifted frames in each direction
        dm = dm + gaussdm(sigma,shift, x-0.5+i*sites, x+0.5+i*sites);
        end
        
    case 'semi-infinite'
        
        for i = -10:11
           dm = dm + gaussdm(sigma,shift, x-0.5+i*sites, x+0.5+i*sites); 
        end
            
        
    otherwise
        error('Error in dispmat: Need to add code to handle different boundary conditions for gaussdisp');
end
        
function d = gaussdm(sigma,shift, neg, pos)


d = normcdf(pos,shift,sigma) - normcdf(neg,shift,sigma);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dm = uniformdisp ( sites, sigma, shift, boundaries )

% generate uniform dispersal matrix

x = repmat(1:sites,[sites,1]);

x = x' - x;
  
dm = zeros( sites );


switch boundaries
    case 'nothing'
        
        dm = unifdm(sigma,shift,x-0.5,x+0.5);
        
    case 'circular'
        
        for i = -10:1:10 % loop over 10 shifted frames in each direction
        dm = dm + gaussdm(sigma,shift, x-0.5+i*sites, x+0.5+i*sites);
        end
        
        
    otherwise 
              error('Error in dispmat: Need to add code to handle different boundary conditions for gaussdisp');
end
  
function d = unifdm(sigma,shift, neg, pos)


d = unifcdf(pos,shift,sigma) - unifcdf(neg,shift,sigma);    

