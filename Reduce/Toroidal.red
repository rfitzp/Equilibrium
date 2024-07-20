% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to find relationship between flux coordinates and toroidal coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linelength 120$

out "Toroidal.out"$

write "<< Output from Toroidal.red >>";

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express inverse of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure inv(ex,var)$

begin scalar c0,c1,c2,i0,i1,i2$
c0 := coeffn(ex,var,0)$
c1 := coeffn(ex,var,1)/c0$
c2 := coeffn(ex,var,2)/c0$
i1 := - c1$
i2 := c1**2 - c2$
return (1 + var*i1 + var**2*i2)/c0

end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express square root of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure sqr(ex,var)$

begin scalar c0,c1,c2,i0,i1,i2$
c0 := coeffn(ex,var,0)$
c1 := coeffn(ex,var,1)/c0$
c2 := coeffn(ex,var,2)/c0$
i1 := c1/2$
i2 := - c1*c1/8 + c2/2$
return c0**(1/2)*(1 + var*i1 + var**2*i2)

end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to average var**2 terms in general polynomial ex
% by performing an integral over general angular variable ivar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure av2(ex,var,ivar)$

begin scalar c0,c1,c2$
c0 := coeffn(ex,var,0)$
c1 := coeffn(ex,var,1)$
c2 := coeffn(ex,var,2)$
c2 := int(c2,ivar)$
return ( sub(ivar=2*pi,c2) - sub(ivar=0,c2) )*var**2/(2*pi)
       + c1*var + c0

end$

% %%%%%%%%%%%%%%%%%%%%%%%
% Useful trig. identities
% %%%%%%%%%%%%%%%%%%%%%%%

for all x, y let cos(x) * cos(y) = ( cos(x+y) + cos(x-y) )/2,
    	     	 cos(x) * sin(y) = ( sin(x+y) - sin(x-y) )/2,
		 sin(x) * sin(y) = ( cos(x-y) - cos(x+y) )/2$
for all x    let sin(x)**2 = ( 1 - cos(2*x) )/2,
                 cos(x)**2 = ( 1 + cos(2*x) )/2$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure that inverse aspect ratio eps is ordered first, both
% internally and externally
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

korder eps,i$
order  eps,i$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure that all quantities expressed as polynomials in eps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor eps, cos, sin$
on  rat$
on  revpri$
off ratpri$

% %%%%%%%%%%%%%%%%%%%%%%%%
% Truncate at order eps**2
% %%%%%%%%%%%%%%%%%%%%%%%%

let eps**3 = 0$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental metric variables:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dh   - lowest-order horizontal shift
% eh   - lowest-order horizontal ellipticity
% th   - lowest-order horizontal triangularity
% sh   - lowest-order horizontal squareness
% dv   - lowest-order vertical shift
% ev   - lowest-order vertical ellipticity
% tv   - lowest-order vertical triangularity
% sv   - lowest-order vertical squareness
% p2   - re-labelling parameter
% x, z - cylindrical coordinates

depend dh,r$
depend eh,r$
depend th,r$
depend sh,r$
depend dv,r$
depend ev,r$
depend tv,r$
depend sv,r$
depend p2,r$
depend x,r,w$
depend z,r,w$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express x and z in terms of flux-surface label r, geometric
% poloidal angle w, and previously defined quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p2 :=   r**3/8
      - r*dh/2
      - (2-1)*eh*eh/(2*r)
      - (3-1)*th*th/(2*r)
      - (4-1)*sh*sh/(2*r)
      - (2-1)*ev*ev/(2*r)
      - (3-1)*tv*tv/(2*r)
      - (4-1)*sv*sv/(2*r)$

x := - r        *cos(w)
     + eps   *dh
     + eps   *eh*cos(w)
     + eps   *th*cos(2*w)
     + eps   *sh*cos(3*w)
     + eps   *ev*sin(w)
     + eps   *tv*sin(2*w)
     + eps   *sv*sin(3*w)
     + eps**2*p2*cos(w)$
z := + r        *sin(w)
     + eps   *eh*sin(w)
     + eps   *th*sin(2*w)
     + eps   *sh*sin(3*w)
     - eps   *dv
     - eps   *ev*cos(w)
     - eps   *tv*cos(2*w)
     - eps   *sv*cos(3*w)
     - eps**2*p2*sin(w)$

% %%%%%%%%%%%%%%%%%%%%%
% Calculate dd1 and dd2
% %%%%%%%%%%%%%%%%%%%%%

x1   := 2 + eps*x$
x2   := x1*x1$
dd1s := (x2 + eps*eps*z*z)/4$
dd1  := sqr(dd1s,eps)$
dd1  := av2(dd1,eps,w)$

%write "d1 := dd10 + eps*dd11 + eps**2*dd12"$
dd10 := coeffn(dd1,eps,0)$
dd11 := coeffn(dd1,eps,1)$
dd12 := coeffn(dd1,eps,2)$

%............................
dd11_target := - r*cos(w)/2$
dd12_target := r*r/16 + dh/2$
%............................

write "dd11 residual:";
dd11 - dd11_target;
write "dd12 residual:";
dd12 - dd12_target;

dd2s := (x*x + z*z)/r/r$

dd2  := sqr(dd2s,eps)$
dd2  := av2(dd2,eps,w)$

%write "d2 := dd20 + eps*dd21 + eps**2*dd22"$
dd20 := coeffn(dd2,eps,0)$
dd21 := coeffn(dd2,eps,1)$
dd22 := coeffn(dd2,eps,2)$

%........................................................................................
dd21_target := - dh*cos(w)/r - eh*cos(2*w)/r - th*cos(3*w)/r - sh*cos(4*w)/r
	       - dv*sin(w)/r - ev*sin(2*w)/r - tv*sin(3*w)/r - sv*sin(4*w)/r$
dd22_target := - r**2/8 + dh/2 + (2*1-1)*dh*dh/4/r/r + (2*2-1)*eh*eh/4/r/r + (2*3-1)*th*th/4/r/r + (2*4-1)*sh*sh/4/r/r
	       	 	       + (2*1-1)*dv*dv/4/r/r + (2*2-1)*ev*ev/4/r/r + (2*3-1)*tv*tv/4/r/r + (2*4-1)*sv*sv/4/r/r$
%........................................................................................

write "dd21 residual:";
dd21 - dd21_target;
write "dd22 residual:";
dd22 - dd22_target;

idd1 := av2( inv(dd1,eps), eps, w)$
idd2 := av2( inv(dd2,eps), eps, w)$

% %%%%%%%%%%%%%
% Calculate eta
% %%%%%%%%%%%%%

coseta := av2( (dd1s - 1 + eps*eps*r*r*dd2s/4)*idd1*idd2/eps/r, eps, w)$

%write "coseta := coseta0 + eps*coseta1";
coseta0 := coeffn(coseta,eps,0)$
coseta1 := coeffn(coseta,eps,1)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express all trig. functions in terms of eta
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%...................................................................................................................
feta := - r*sin(eta)/2 
     	+ cos(pi)*dh*sin(eta)/r + cos(2*pi)*eh*sin(2*eta)/r + cos(3*pi)*th*sin(3*eta)/r + cos(4*pi)*sh*sin(4*eta)/r
        + cos(pi)*dv*cos(eta)/r + cos(2*pi)*ev*cos(2*eta)/r + cos(3*pi)*tv*cos(3*eta)/r + cos(4*pi)*sv*cos(4*eta)/r$
%...................................................................................................................	

let           cos(w)   =           - (  cos(eta)   + sin(eta)  *eps  *feta)$
for all n let cos(n*w) = cos(n*pi) * (  cos(n*eta) + sin(n*eta)*eps*n*feta)$
let           sin(w)   =           - (- sin(eta)   + cos(eta)  *eps  *feta)$
for all n let sin(n*w) = cos(n*pi) * (- sin(n*eta) + cos(n*eta)*eps*n*feta)$

coseta0 := coeffn(coseta,eps,0);
coseta1 := coeffn(coseta,eps,1);

% %%%%%%%%%%%%
% Calculate xi
% %%%%%%%%%%%%

xi := av2( (dd1s + eps*eps*r*r*dd2s/4)*idd1*idd2, eps, eta)$

%write "xi := xi0 + eps*xi1 + eps**2*xi2";
xi0 := coeffn(xi,eps,0)$
xi1 := coeffn(xi,eps,1)$
xi2 := coeffn(xi,eps,2)$

%...........................................................................................
xi1_target := r*cos(eta)/2 + cos(pi)*dh*cos(eta)/r + cos(2*pi)*eh*cos(2*eta)/r + cos(3*pi)*th*cos(3*eta)/r + cos(4*pi)*sh*cos(4*eta)/r
	                   - cos(pi)*dv*sin(eta)/r - cos(2*pi)*ev*sin(2*eta)/r - cos(3*pi)*tv*sin(3*eta)/r - cos(4*pi)*sv*sin(4*eta)/r$	  
xi2_target := 5*r*r/16 - dh/4 + 3*dh*dh/4/r/r + 3*eh*eh/4/r/r + 3*th*th/4/r/r + 3*sh*sh/4/r/r
 	      	       	      + 3*dv*dv/4/r/r + 3*ev*ev/4/r/r + 3*tv*tv/4/r/r + 3*sv*sv/4/r/r$
%............................................................................................

write "xi1 residual:";
xi1 - xi1_target;
write "xi2 residual:";
xi2 - xi2_target;

% %%%%%%%%%%%%%%%%%
% Calculate log(xi)
% %%%%%%%%%%%%%%%%%

xired := eps*xi1 + eps**2*xi2$
lnxi  := av2( xired - xired*xired/2, eps, eta)$

%write "ln(xi) := lnxi0 + eps*lnxi1 + eps**2*lnxi2";
lnxi0 := coeffn(lnxi,eps,0)$
lnxi1 := coeffn(lnxi,eps,1)$
lnxi2 := coeffn(lnxi,eps,2)$

%.............................................................................................
lnxi1_target := r*cos(eta)/2 + cos(pi)*dh*cos(eta)/r + cos(2*pi)*eh*cos(2*eta)/r + cos(3*pi)*th*cos(3*eta)/r + cos(4*pi)*sh*cos(4*eta)/r
	     		     - cos(pi)*dv*sin(eta)/r - cos(2*pi)*ev*sin(2*eta)/r - cos(3*pi)*tv*sin(3*eta)/r - cos(4*pi)*sv*sin(4*eta)/r$
lnxi2_target := r*r/4 + dh*dh/2/r/r + eh*eh/2/r/r + th*th/2/r/r + sh*sh/2/r/r
	     	      + dv*dv/2/r/r + ev*ev/2/r/r + tv*tv/2/r/r + sv*sv/2/r/r$
%.............................................................................................

write "lnxi1 residual:";
lnxi1 - lnxi1_target;
write "lnxi2 residual:";
lnxi2 - lnxi2_target;

out t$

;bye;
