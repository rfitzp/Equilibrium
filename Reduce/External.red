% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to find relationship between cylindrical coordinates and toroidal coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linelength 120$

out "External.out"$

write "<< Output from External.red >>";

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express nth power of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure pow1(ex,var,n)$

begin scalar c0,c1,c2,i0,i1,i2$
c0 := coeffn(ex,var,0)$
c1 := coeffn(ex,var,1)/c0$
c2 := coeffn(ex,var,2)/c0$
i1 := n*c1$
i2 := n*(n-1)*c1*c1/2 + n*c2$
return c0**n*(1 + var*i1 + var**2*i2)

end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to express n/m th power of 2nd order polynomial ex (in var)
% as another 2nd order polynomial (in var)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure pow2(ex,var,n,m)$

begin scalar c0,c1,c2,i0,i1,i2$
c0 := coeffn(ex,var,0)$
c1 := coeffn(ex,var,1)/c0$
c2 := coeffn(ex,var,2)/c0$
i1 := n*c1/m$
i2 := n*(n-m)*c1*c1/2/m/m + n*c2/m$
return c0**(n/m)*(1 + var*i1 + var**2*i2)

end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to return amplitude of nth-order Fourier cosine
% component in angle-like variable ivar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure fcos(ex,ivar,n)$

begin scalar x,y,norm$

x := ex*cos(n*ivar)$
y := int(x,ivar)$
if n=0 then norm := 1 else norm := 2$
return norm*( sub(ivar=2*pi,y) - sub(ivar=0,y) )/(2*pi)

end$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to return amplitude of nth-order Fourier sine
% component in angle-like variable ivar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

procedure fsin(ex,ivar,n)$

begin scalar x,y,norm$

x := ex*sin(n*ivar)$
y := int(x,ivar)$
if n=0 then norm := 1 else norm := 2$
return norm*( sub(ivar=2*pi,y) - sub(ivar=0,y) )/(2*pi)

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

depend x,r,w$
depend z,r,w$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express x and z in terms of radius r, and geometric
% poloidal angle w
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x := r * cos(w)$
z := r * sin(w)$

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

%.........................
dd11_target := r*cos(w)/2$
dd12_target := r*r/16$
%.........................

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

%................
dd20_target := 1$
dd21_target := 0$
dd22_target := 0$
%................

write "dd20 residual:";
dd20 - dd20_target;
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

%....................
feta := r*sin(eta)/2$
%....................

let           cos(w)   = (cos(eta)   - sin(eta)  *eps  *feta)$
for all n let cos(n*w) = (cos(n*eta) - sin(n*eta)*eps*n*feta)$
let           sin(w)   = (sin(eta)   + cos(eta)  *eps  *feta)$
for all n let sin(n*w) = (sin(n*eta) + cos(n*eta)*eps*n*feta)$

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

%..........................
xi0_target := 0$
xi1_target := r*cos(eta)/2$	  
xi2_target := 3*r*r/16$
%..........................

write "xi0 residual:";
xi0 - xi0_target;
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

%............................
lnxi0_target := 0$
lnxi1_target := r*cos(eta)/2$  
lnxi2_target := r*r/8$
%............................

write "lnxi0 residual:";
lnxi0 - lnxi0_target;
write "lnxi1 residual:";
lnxi1 - lnxi1_target;
write "lnxi2 residual:";
lnxi2 - lnxi2_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Poloidal flux from far-vacuum region
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi0 := 1$
xi1 := r*cos(eta)/2$
xi2 := 3*r*r/16$

xi := xi0 + eps*xi1 + eps**2*xi2$

xir := xi/r$
rxi := r * av2( inv(xi,eps), eps, eta)$

fac     := 1 - eps*rxi*cos(eta)$
factor1 := av2( pow2(fac,eps,-1,2), eps, eta)$

lnxi0 := 0$
lnxi1 := r*cos(eta)/2$
lnxi2 := r*r/8$

lnxi := lnxi0 + eps*lnxi1 + eps**2*lnxi2$

factor2 :=   qc00 * (1 - (eps***2)*(rxi**2)/16)
	   + qc1 * eps*rxi     *cos(eta)
	   + qc2 * eps*(rxi**2)*cos(2*eta)
	   + qc3 * eps*(rxi**3)*cos(3*eta)
	   + qc4 * eps*(rxi**4)*cos(4*eta)
	   + qs1 * eps*rxi     *sin(eta)
	   + qs2 * eps*(rxi**2)*sin(2*eta)
	   + qs3 * eps*(rxi**3)*sin(3*eta)
	   + qs4 * eps*(rxi**4)*sin(4*eta)
	   + qc02 * eps*eps$

factor2 = av2(factor2,eps,eta)$

psivac := av2(factor1*factor2,eps,eta)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express all trig. functions in terms of t
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%................
ft := r*sin(t)/2$
%................

let           cos(eta)   = (cos(t)   + sin(t)  *eps  *ft)$
for all n let cos(n*eta) = (cos(n*t) + sin(n*t)*eps*n*ft)$
let           sin(eta)   = (sin(t)   - cos(t)  *eps  *ft)$
for all n let sin(n*eta) = (sin(n*t) - cos(n*t)*eps*n*ft)$

psivac := av2(factor1*factor2,eps,t)$

psivac0   := coeffn(psivac,eps,0)$
psivac1   := coeffn(psivac,eps,1)$
psivac1c0 := fcos(psivac1,t,0)$
psivac1c1 := fcos(psivac1,t,1)$
psivac1c2 := fcos(psivac1,t,2)$
psivac1c3 := fcos(psivac1,t,3)$
psivac1c4 := fcos(psivac1,t,4)$
psivac1s1 := fsin(psivac1,t,1)$
psivac1s2 := fsin(psivac1,t,2)$
psivac1s3 := fsin(psivac1,t,3)$
psivac1s4 := fsin(psivac1,t,4)$
psivac2   := coeffn(psivac,eps,2)$

%....................................
psivac0_target   := qc00$

psivac1c1_target := (qc00/2 + qc1)*r$
psivac1c2_target := qc2*r*r$
psivac1c3_target := qc3*r*r*r$
psivac1c4_target := qc4*r*r*r*r$

psivac1s1_target := qs1*r$
psivac1s2_target := qs2*r*r$
psivac1s3_target := qs3*r*r*r$
psivac1s4_target := qs4*r*r*r*r$

psivac2_target   := (qc00/2 + qc1)*r*r/4 + qc02$
%...................................
		  
write "psivac0 residual:";
psivac0 - psivac0_target;

write "psivac1c1 residual:";
psivac1c1 - psivac1c1_target;
write "psivac1c2 residual:";
psivac1c2 - psivac1c2_target;
write "psivac1c3 residual:";
psivac1c3 - psivac1c3_target;
write "psivac1c4 residual:";
psivac1c4 - psivac1c4_target;

write "psivac1s1 residual:";
psivac1s1 - psivac1s1_target;
write "psivac1s2 residual:";
psivac1s2 - psivac1s2_target;
write "psivac1s3 residual:";
psivac1s3 - psivac1s3_target;
write "psivac1s4 residual:";
psivac1s4 - psivac1s4_target;

write "psivac2 residual:";
psivac2 - psivac2_target;

out t$

;bye;
