% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to characterize aspect-ratio-expanded tokamak equilibrium
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linelength 120$

out "Metric.out"$

write "<< Output from Metric.red >>";

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

% dh   - lowest-order horizontal shift
% eh   - lowest-order horizontal ellipticity
% th   - lowest-order horizontal triangularity
% sh   - lowest-order horizontal squareness
% dv   - lowest-order vertical shift
% ev   - lowest-order vertical ellipticity
% tv   - lowest-order vertical triangularity
% sv   - lowest-order vertical squareness
% l2   - re-labelling parameter
% x, z - cylindrical coordinates

depend dh,r$
depend eh,r$
depend th,r$
depend sh,r$
depend dv,r$
depend ev,r$
depend tv,r$
depend sv,r$
depend l2,r$
depend x,r,w$
depend z,r,w$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express x and z in terms of flux-surface label r, geometric
% poloidal angle w, and previously defined quantities
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x := - r        *cos(w)
     + eps   *dh
     + eps   *eh*cos(w)
     + eps   *th*cos(2*w)
     + eps   *sh*cos(3*w)
     + eps   *ev*sin(w)
     + eps   *tv*sin(2*w)
     + eps   *sv*sin(3*w)
     + eps**2*l2*cos(w)$
z := + r        *sin(w)
     + eps   *eh*sin(w)
     + eps   *th*sin(2*w)
     + eps   *sh*sin(3*w)
     - eps   *dv
     - eps   *ev*cos(w)
     - eps   *tv*cos(2*w)
     - eps   *sv*cos(3*w)
     - eps**2*l2*sin(w)$

%write "R := xx0 + eps*xx1 + eps**2*xx2 + eps**3*xx3"$
xx0 := 1$
xx1 := coeffn(x,eps,0)$
xx2 := coeffn(x,eps,1)$
xx3 := coeffn(x,eps,2)$

%write "Z := zz0 + eps*zz1 + eps**2*zz2 + eps**3*zz3"$
zz0 := 0$
zz1 := coeffn(z,eps,0)$
zz2 := coeffn(z,eps,1)$
zz3 := coeffn(z,eps,2)$

% %%%%%%%%%%%%%%%%%
% Evaluate Jacobian
% %%%%%%%%%%%%%%%%%

xr := df(x,r)$
xw := df(x,w)$
zr := df(z,r)$
zw := df(z,w)$
j  := xw*zr - xr*zw$

%write "J := j0 + eps*j1 + eps**2*j2"$
j0 := coeffn(j,eps,0)$
j1 := coeffn(j,eps,1)$
j2 := coeffn(j,eps,2)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform to flux coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1  := 1 + eps*x$
ix1 := inv(x1,eps)$
k   := j*ix1$
ik  := int(k,w)$
ika := ( sub(w=2*pi,ik) - sub(w=0,ik) ) /(2*pi)$

%write "J/R := k0 + eps*k1 + eps**2*k2"$
k0 := coeffn(k,eps,0)$
k1 := coeffn(k,eps,1)$
k2 := coeffn(k,eps,2)$

%............................
l2_target :=   r**3/8
             - r*dh/2
	     - (2-1)*eh*eh/(2*r)
	     - (3-1)*th*th/(2*r)
	     - (4-1)*sh*sh/(2*r)
	     - (2-1)*ev*ev/(2*r)
	     - (3-1)*tv*tv/(2*r)
	     - (4-1)*sv*sv/(2*r)$
%............................

%write "<J/R> := ik0 + eps*ik1 + eps*ik2"%
l2  := l2_target$
ik0 := coeffn(ika,eps,0)$
ik1 := coeffn(ika,eps,1)$
ik2 := coeffn(ika,eps,2)$

%...............
ik0_target := r$
ik1_target := 0$
ik2_target := 0$
%...............

write "ik0 residual:";
ik0 - ik0_target;
write "ik1 residual:";
ik1 - ik1_target;
write "ik2 residual:";
ik2 - ik2_target;

%write "theta := tt0 + eps*tt1 + eps**2*tt2"$
theta := av2(ik/r, eps, w)$
tt0   := coeffn(theta,eps,0)$
tt1   := coeffn(theta,eps,1)$
tt2   := coeffn(theta,eps,2)$

%..................................................
tt1_target := r*sin(w)
 	      - (df(dh,r) - (1-1)*dh/r)*sin(w)
              - (df(eh,r) - (2-1)*eh/r)*sin(2*w)/2
	      - (df(th,r) - (3-1)*th/r)*sin(3*w)/3
	      - (df(sh,r) - (4-1)*sh/r)*sin(4*w)/4
	      + (df(dv,r) - (1-1)*dv/r)*cos(w)
	      + (df(ev,r) - (2-1)*ev/r)*cos(2*w)/2
	      + (df(tv,r) - (3-1)*tv/r)*cos(3*w)/3
	      + (df(sv,r) - (4-1)*sv/r)*cos(4*w)/4$
%..................................................

write "tt1 residual:";
tt1 - tt1_target;

write "Secular variation residual:";
( sub(w=2*pi,theta) - sub(w=0,theta) )/(2*pi) - 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate covariant metric tensor elements in geometric coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gww := zw**2 + xw**2$
grw := zr*zw + xr*xw$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate contravariant metric tensor elements in geometric coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ij     := inv(j,eps)$
ij2    := ij*ij$
gradr2 :=   gww*ij2$
grrgrw := - grw*ij2$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate contravariant metric tensor elements in flux coordinates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tw     := df(theta,w)$
ttr    := df(theta,r)$
grrgrt := grrgrw*tw + gradr2*ttr$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate geometric angle w as function of flux-coordinate angle t
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only need transformation to order eps.

% t = w + eps*ww1(w) + O(eps**2)

% t = w + eps*ww1(t) + O(eps**2)

% w = t + eps*tth(t) + O(eps**2)

% where tth(w) = - ww1(t)

www1 := coeffn(theta,eps,1)$
ttth := - sub(w=t,www1)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express all trig. functions in terms of t
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

let           cos(w)   = cos(t)   - sin(t)*eps*ttth$
for all n let cos(n*w) = cos(n*t) - n*sin(n*t)*eps*ttth$
let           sin(w)   = sin(t)   + cos(t)*eps*ttth$
for all n let sin(n*w) = sin(n*t) + n*cos(n*t)*eps*ttth$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output metric tensor elements with t-averaged eps**2 coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%write "gr2 := gr20 + eps*gr21 + eps**2*gr22"$
gr2  := av2(gradr2,eps,t)$
gr20 := coeffn(gr2,eps,0)$
gr21 := coeffn(gr2,eps,1)$
gr22 := coeffn(gr2,eps,2)$

%.............................................................................
gr21_target := + 2*df(dh,r)*cos(t)
               + 2*df(eh,r)*cos(2*t)
	       + 2*df(th,r)*cos(3*t)
	       + 2*df(sh,r)*cos(4*t)
	       + 2*df(dv,r)*sin(t)
               + 2*df(ev,r)*sin(2*t)
	       + 2*df(tv,r)*sin(3*t)
	       + 2*df(sv,r)*sin(4*t)$
gr22_target := 3*r**2/4
	       - dh
	       + (df(dh,r)*df(dh,r) + (1*1-1)*dh*dh/r/r)/2
	       + (df(eh,r)*df(eh,r) + (2*2-1)*eh*eh/r/r)/2
	       + (df(th,r)*df(th,r) + (3*3-1)*th*th/r/r)/2
	       + (df(sh,r)*df(sh,r) + (4*4-1)*sh*sh/r/r)/2
	       + (df(dv,r)*df(dv,r) + (1*1-1)*dv*dv/r/r)/2
	       + (df(ev,r)*df(ev,r) + (2*2-1)*ev*ev/r/r)/2
	       + (df(tv,r)*df(tv,r) + (3*3-1)*tv*tv/r/r)/2
	       + (df(sv,r)*df(sv,r) + (4*4-1)*sv*sv/r/r)/2$
%.............................................................................

write "gr21 residual:";
gr21 - gr21_target;
write "gr22 residual:";
gr22 - gr22_target;

%write "grt2 := grt20 + eps*grt21 + eps**2*grt22"$
grt2  := av2(grrgrt,eps,t)$
grt20 := coeffn(grt2,eps,0)$
grt21 := coeffn(grt2,eps,1)$
grt22 := coeffn(grt2,eps,2)$

%.......................................................................
grt21_target := sin(t)  
                - (df(dh,r,2) + df(dh,r)/r + (1*1-1)*dh/r/r)*sin(t)  /1             
	     	- (df(eh,r,2) + df(eh,r)/r + (2*2-1)*eh/r/r)*sin(2*t)/2
		- (df(th,r,2) + df(th,r)/r + (3*3-1)*th/r/r)*sin(3*t)/3
		- (df(sh,r,2) + df(sh,r)/r + (4*4-1)*sh/r/r)*sin(4*t)/4
		+ (df(dv,r,2) + df(dv,r)/r + (1*1-1)*dv/r/r)*cos(t)  /1
	     	+ (df(ev,r,2) + df(ev,r)/r + (2*2-1)*ev/r/r)*cos(2*t)/2
		+ (df(tv,r,2) + df(tv,r)/r + (3*3-1)*tv/r/r)*cos(3*t)/3
		+ (df(sv,r,2) + df(sv,r)/r + (4*4-1)*sv/r/r)*cos(4*t)/4$
grt22_target :=   3*df(dv,r)/2 + df(dv,r,2)*r/2
                + (2-1)*df(ev,r)*eh/2/r/r       - (2-1)*df(eh,r)*ev/2/r/r
		+ (3-1)*df(tv,r)*th/2/r/r       - (3-1)*df(th,r)*tv/2/r/r
		+ (4-1)*df(sv,r)*sh/2/r/r       - (4-1)*df(sh,r)*sv/2/r/r
		+ (2-1)*df(ev,r,2)*eh/2/2/r     - (2-1)*df(eh,r,2)*ev/2/2/r 
		+ (3-1)*df(tv,r,2)*th/2/3/r     - (3-1)*df(th,r,2)*tv/2/3/r
		+ (4-1)*df(sv,r,2)*sh/2/4/r     - (4-1)*df(sh,r,2)*sv/2/4/r
		+       df(dv,r,2)*df(dh,r)/2/1 - df(dh,r,2)*df(dv,r)/2/1
		+       df(ev,r,2)*df(eh,r)/2/2 - df(eh,r,2)*df(ev,r)/2/2 
		+       df(tv,r,2)*df(th,r)/2/3 - df(th,r,2)*df(tv,r)/2/3
		+       df(sv,r,2)*df(sh,r)/2/4 - df(sh,r,2)*df(sv,r)/2/4$
%.......................................................................

write "grt21 residual:";
grt21 - grt21_target;
write "grt22 residual:";
grt22 - grt22_target;

%write "x*x := x20 + eps*x21 + eps**2*x22"$
x1  := 1 + eps*x$
x2  := x1*x1$
x2  := av2(x2,eps,t)$
x20 := coeffn(x2,eps,0)$
x21 := coeffn(x2,eps,1)$
x22 := coeffn(x2,eps,2)$

%...........................................
x21_target := - 2*r*cos(t)$
x22_target := - (r*r/2 - 2*dh - r*df(dh,r))$
%...........................................

write "x21 residual:";
x21 - x21_target;
write "x22 residual:";
x22 - x22_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental Grad-Shafranov equation variables:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f1 - lowest-order poloidal flux function
% f3 - higher-order poloidal flux function
% g2 - lowest-order toroidal flux function
% g4 - higher-order toroidal flux function
% p2 - lowest-order pressure
% p4 - higher-order pressure

depend f1,r$
depend f3,r$
depend g2,r$
depend g4,r$
depend p2,r$
depend p4,r$

ff  := f1 + eps**2*f3$
gg  := 1 + eps**2*g2$
ggp := df(g2,r) + eps*eps*df(g4,r)$
pp  := df(p2,r) + eps*eps*df(p4,r)$

% %%%%%%%%%%%%%%%%%%%%%%%
% Grad-Shafranov equation
% %%%%%%%%%%%%%%%%%%%%%%%

fac1 := ff * gr2$
fac2 := ff * grt2$

gs := ff * ( df(fac1,r) + df(fac2,t) ) + r*r * (gg*ggp + x2*pp)$
gs := av2(gs,eps,t)$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier analyze Grad-Shafranov equation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gs0  := fcos(gs,t,0)$
gs00 := coeffn( gs0, eps, 0)$
gs02 := coeffn( gs0, eps, 2)$
gs1c := coeffn( fcos(gs,t,1), eps, 1)$
gs2c := coeffn( fcos(gs,t,2), eps, 1)$
gs3c := coeffn( fcos(gs,t,3), eps, 1)$
gs4c := coeffn( fcos(gs,t,4), eps, 1)$
gs1s := coeffn( fsin(gs,t,1), eps, 1)$
gs2s := coeffn( fsin(gs,t,2), eps, 1)$
gs3s := coeffn( fsin(gs,t,3), eps, 1)$
gs4s := coeffn( fsin(gs,t,4), eps, 1)$

%write "Zeroth-order equilibrium equation:";
c0  := coeffn( gs00, df(g2,r), 0)$
c1  := coeffn( gs00, df(g2,r), 1)$
g2p := -c0/c1$

%write "Shafranov shift equation:";
c0   := coeffn( gs1c, df(dh,r,2), 0)$
c1   := coeffn( gs1c, df(dh,r,2), 1)$
dhpp := -c0/c1$

%write "Horizontal ellpiticity equation:";
c0   := coeffn( gs2c, df(eh,r,2), 0)$
c1   := coeffn( gs2c, df(eh,r,2), 1)$
ehpp := -c0/c1$

%write "Horizontal triangularity equation:";
c0   := coeffn( gs3c, df(th,r,2), 0)$
c1   := coeffn( gs3c, df(th,r,2), 1)$
thpp := -c0/c1$

%write "Horizontal squareness equation:";
c0   := coeffn( gs4c, df(sh,r,2), 0)$
c1   := coeffn( gs4c, df(sh,r,2), 1)$
shpp := -c0/c1$

%write "Vertical shift equation:";
c0   := coeffn( gs1s, df(dv,r,2), 0)$
c1   := coeffn( gs1s, df(dv,r,2), 1)$
dvpp := -c0/c1$

%write "Vertical ellpiticity equation:";
c0   := coeffn( gs2s, df(ev,r,2), 0)$
c1   := coeffn( gs2s, df(ev,r,2), 1)$
evpp := -c0/c1$

%write "Vertical triangularity equation:";
c0   := coeffn( gs3s, df(tv,r,2), 0)$
c1   := coeffn( gs3s, df(tv,r,2), 1)$
tvpp := -c0/c1$

%write "Vertical squareness equation:";
c0   := coeffn( gs4s, df(sv,r,2), 0)$
c1   := coeffn( gs4s, df(sv,r,2), 1)$
svpp := -c0/c1$

%write "Second-order equilibrium equation:";
c0  := coeffn( gs02, df(f3,r), 0)$
c1  := coeffn( gs02, df(f3,r), 1)$
f3p := -c0/c1$
f3p := sub(df(g2,r)  =g2p, f3p)$
f3p := sub(df(dh,r,2)=dhpp,f3p)$
f3p := sub(df(eh,r,2)=ehpp,f3p)$
f3p := sub(df(th,r,2)=thpp,f3p)$
f3p := sub(df(sh,r,2)=shpp,f3p)$
f3p := sub(df(dv,r,2)=dvpp,f3p)$
f3p := sub(df(ev,r,2)=evpp,f3p)$
f3p := sub(df(tv,r,2)=tvpp,f3p)$
f3p := sub(df(sv,r,2)=svpp,f3p)$

ff3 := - f1*f1 * (3*r*r/2 - 2*r*df(dh,r)
               + (df(dh,r)*df(dh,r) + 2*(1*1-1)*df(dh,r)*eh/r - (1*1-1)*dh*dh/r/r)
       	       + (df(eh,r)*df(eh,r) + 2*(2*2-1)*df(eh,r)*eh/r - (2*2-1)*eh*eh/r/r)
	       + (df(th,r)*df(th,r) + 2*(3*3-1)*df(th,r)*th/r - (3*3-1)*th*th/r/r)
	       + (df(sh,r)*df(sh,r) + 2*(4*4-1)*df(sh,r)*sh/r - (4*4-1)*sh*sh/r/r)
	       + (df(dv,r)*df(dv,r) + 2*(1*1-1)*df(dv,r)*dv/r - (1*1-1)*dv*dv/r/r)
	       + (df(ev,r)*df(ev,r) + 2*(2*2-1)*df(ev,r)*ev/r - (2*2-1)*ev*ev/r/r)
	       + (df(tv,r)*df(tv,r) + 2*(3*3-1)*df(tv,r)*tv/r - (3*3-1)*tv*tv/r/r)
	       + (df(sv,r)*df(sv,r) + 2*(4*4-1)*df(sv,r)*sv/r - (4*4-1)*sv*sv/r/r))/r
	+ df(f1,r)*f1 * (g2 - 3*r*r/4 + dh
	                + (3*df(dh,r)*df(dh,r) - (1*1-1)*dh*dh/r/r)/2
	  	      	+ (3*df(eh,r)*df(eh,r) - (2*2-1)*eh*eh/r/r)/2
			+ (3*df(th,r)*df(th,r) - (3*3-1)*th*th/r/r)/2
			+ (3*df(sh,r)*df(sh,r) - (4*4-1)*sh*sh/r/r)/2
			+ (3*df(dv,r)*df(dv,r) - (1*1-1)*dv*dv/r/r)/2
			+ (3*df(ev,r)*df(ev,r) - (2*2-1)*ev*ev/r/r)/2
			+ (3*df(tv,r)*df(tv,r) - (3*3-1)*tv*tv/r/r)/2
			+ (3*df(sv,r)*df(sv,r) - (4*4-1)*sv*sv/r/r)/2)
	+ r*r*df(p2,r) * (g2 + r*r/2 - 2*dh - 3*r*df(dh,r))$		

%...........................................................................
gs_target   := gs00
	       + eps * (gs1c*cos(t) + gs2c*cos(2*t) + gs3c*cos(3*t) + gs4c*cos(4*t))
	       + eps * (gs1s*sin(t) + gs2s*sin(2*t) + gs3s*sin(3*t) + gs4s*sin(4*t))
 	       + eps*eps * gs02$
g2p_target  := - df(p2,r) - f1 * df(f1,r)/r**2$
dhpp_target := - (2*df(f1,r)/f1 - 1/r)*df(dh,r) - 1 + 2*r**3*df(p2,r)/f1**2$
ehpp_target := - (2*df(f1,r)/f1 - 1/r)*df(eh,r) + (2*2-1)*eh/r**2$
thpp_target := - (2*df(f1,r)/f1 - 1/r)*df(th,r) + (3*3-1)*th/r**2$
shpp_target := - (2*df(f1,r)/f1 - 1/r)*df(sh,r) + (4*4-1)*sh/r**2$
dvpp_target := - (2*df(f1,r)/f1 - 1/r)*df(dv,r)$
evpp_target := - (2*df(f1,r)/f1 - 1/r)*df(ev,r) + (2*2-1)*ev/r**2$
tvpp_target := - (2*df(f1,r)/f1 - 1/r)*df(tv,r) + (3*3-1)*tv/r**2$
svpp_target := - (2*df(f1,r)/f1 - 1/r)*df(sv,r) + (4*4-1)*sv/r**2$
f3p_target  := - df(f1,r)*f3/f1 + ff3/f1$
%...........................................................................

write "g2p residual:";
g2p - g2p_target;
write "dhpp residual:";
dhpp - dhpp_target;
write "ehpp residual:";
ehpp - ehpp_target;
write "thpp residual:";
thpp - thpp_target;
write "shpp residual:";
shpp - shpp_target;
write "dvpp residual:";
dvpp - dvpp_target;
write "evpp residual:";
evpp - evpp_target;
write "tvpp residual:";
tvpp - tvpp_target;
write "svpp residual:";
svpp - svpp_target;
write "f3p residual:";
f3p - f3p_target;
write "gs residual:";
gs - gs_target;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental vacuum solution variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f1a  - value of f1  at plasma boundary
% f3a  - value of f3  at plasma boundary
% dha  - value of dh  at plasma boundary
% dhpa - value of dh' at plasma boundary
% eha  - value of eh  at plasma boundary
% ehpa - value of eh' at plasma boundary
% tha  - value of th  at plasma boundary
% thpa - value of th' at plasma boundary
% sha  - value of sh  at plasma boundary
% shpa - value of sh' at plasma boundary
% dva  - value of dv  at plasma boundary
% dvpa - value of dv' at plasma boundary
% eva  - value of ev  at plasma boundary
% evpa - value of ev' at plasma boundary
% tva  - value of tv  at plasma boundary
% tvpa - value of tv' at plasma boundary
% sva  - value of sv  at plasma boundary
% svpa - value of sv' at plasma boundary

% f1v - f1 in vacuum
% dhv - dh in vacuum
% ehv - eh in vacuum
% thv - th in vacuum
% shv - sh in vacuum
% dvv - dv in vacuum
% evv - ev in vacuum
% tvv - tv in vacuum
% svv - sv in vacuum

% psi1 - lowest-order poloidal flux
% psi3 - higher-order poloidal flux

% psi1a - value of psi1 on plasma boundary
% psi3a - value of psi3 on plasma boundary

depend dhv,r$
depend ehv,r$
depend thv,r$
depend shv,r$
depend dvv,r$
depend evv,r$
depend tvv,r$
depend svv,r$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate vacuum poloidal flux
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1v := f1a$
g2v := 0$

ff3v := - f1v*f1v * (3*r*r/2 - 2*r*df(dhv,r)
                  + (df(dhv,r)*df(dhv,r) + 2*(1*1-1)*df(dhv,r)*dhv/r - (1*1-1)*dhv*dhv/r/r)
       	     	  + (df(ehv,r)*df(ehv,r) + 2*(2*2-1)*df(ehv,r)*ehv/r - (2*2-1)*ehv*ehv/r/r)
	       	  + (df(thv,r)*df(thv,r) + 2*(3*3-1)*df(thv,r)*thv/r - (3*3-1)*thv*thv/r/r)
	       	  + (df(shv,r)*df(shv,r) + 2*(4*4-1)*df(shv,r)*shv/r - (4*4-1)*shv*shv/r/r)
	     	  + (df(dvv,r)*df(dvv,r) + 2*(1*1-1)*df(dvv,r)*dvv/r - (1*1-1)*dvv*dvv/r/r)
	          + (df(evv,r)*df(evv,r) + 2*(2*2-1)*df(evv,r)*evv/r - (2*2-1)*evv*evv/r/r)
	       	  + (df(tvv,r)*df(tvv,r) + 2*(3*3-1)*df(tvv,r)*tvv/r - (3*3-1)*tvv*tvv/r/r)
	       	  + (df(svv,r)*df(svv,r) + 2*(4*4-1)*df(svv,r)*svv/r - (4*4-1)*svv*svv/r/r))/r
	+ df(f1v,r)*f1v * (g2v - 3*r*r/4 + dhv
	                + (3*df(dhv,r)*df(dhv,r) - (1*1-1)*dhv*dhv/r/r)/2
	  	      	+ (3*df(ehv,r)*df(ehv,r) - (2*2-1)*ehv*ehv/r/r)/2
			+ (3*df(thv,r)*df(thv,r) - (3*3-1)*thv*thv/r/r)/2
			+ (3*df(shv,r)*df(shv,r) - (4*4-1)*shv*shv/r/r)/2
			+ (3*df(dvv,r)*df(dvv,r) - (1*1-1)*dvv*dvv/r/r)/2
			+ (3*df(evv,r)*df(evv,r) - (2*2-1)*evv*evv/r/r)/2
			+ (3*df(tvv,r)*df(tvv,r) - (3*3-1)*tvv*tvv/r/r)/2
			+ (3*df(svv,r)*df(svv,r) - (4*4-1)*svv*svv/r/r)/2)
	+ r*r*df(p2v,r) * (g2v + r*r/2 - 2*dhv - 3*r*df(dhv,r))$		


ff3pv := - df(f1v,r)*f3v/f1v + ff3v/f1v$

%..................................................................................................
ff3pv_target := - f1a * (3*r*r/2 - 2*r*df(dhv,r)
	     	      	+ df(dhv,r)*df(dhv,r) + 2*(1*1-1)*df(dhv,r)*dhv/r - (1*1-1)*dhv*dhv/r/r
	    	     	+ df(ehv,r)*df(ehv,r) + 2*(2*2-1)*df(ehv,r)*ehv/r - (2*2-1)*ehv*ehv/r/r
			+ df(thv,r)*df(thv,r) + 2*(3*3-1)*df(thv,r)*thv/r - (3*3-1)*thv*thv/r/r
			+ df(shv,r)*df(shv,r) + 2*(4*4-1)*df(shv,r)*shv/r - (4*4-1)*shv*shv/r/r
			+ df(dvv,r)*df(dvv,r) + 2*(1*1-1)*df(dvv,r)*dvv/r - (1*1-1)*dvv*dvv/r/r
	    	     	+ df(evv,r)*df(evv,r) + 2*(2*2-1)*df(evv,r)*evv/r - (2*2-1)*evv*evv/r/r
			+ df(tvv,r)*df(tvv,r) + 2*(3*3-1)*df(tvv,r)*tvv/r - (3*3-1)*tvv*tvv/r/r
			+ df(svv,r)*df(svv,r) + 2*(4*4-1)*df(svv,r)*svv/r - (4*4-1)*svv*svv/r/r)/r$
%..................................................................................................

write "ff3pv residual:";
ff3pv - ff3pv_target;

ehp := (ehpa - (2+1)*eha)/(2*2)$
ehm := (ehpa + (2-1)*eha)/(2*2)$
thp := (thpa - (3+1)*tha)/(2*3)$
thm := (thpa + (3-1)*tha)/(2*3)$
shp := (shpa - (4+1)*sha)/(2*4)$
shm := (shpa + (4-1)*sha)/(2*4)$
dvp := (dvpa - (1+1)*dva)/(2*1)$
dvm := (dvpa + (1-1)*dva)/(2*1)$
evp := (evpa - (2+1)*eva)/(2*2)$
evm := (evpa + (2-1)*eva)/(2*2)$
tvp := (tvpa - (3+1)*tva)/(2*3)$
tvm := (tvpa + (3-1)*tva)/(2*3)$
svp := (svpa - (4+1)*sva)/(2*4)$
svm := (svpa + (4-1)*sva)/(2*4)$

dhv := dha - r*r*log(r)/2 + (2*dhpa + 1)*(r*r - 1)/4$
ehv := ehm * r**(2+1) - ehp * r**(1-2)$
thv := thm * r**(3+1) - thp * r**(1-3)$
shv := shm * r**(4+1) - shp * r**(1-4)$
dvv := dvm * r**(1+1) - dvp * r**(1-1)$
evv := evm * r**(2+1) - evp * r**(1-2)$
tvv := tvm * r**(3+1) - tvp * r**(1-3)$
svv := svm * r**(4+1) - svp * r**(1-4)$

ff3v := int(ff3pv,r)$
ff3v := f3a + ff3v - sub(r=1,ff3v)$

%..............................................................................................
ff3v_target :=   f3a
	       - f1a * ((1-2*dhpa)*r*r*log(r) + r*r*log(r)*log(r) + (1-dhpa+dhpa*dhpa)*(r*r-1)
	     	     + (2*2)*(2+1)*ehm*ehm*(r**(2*2)-1) + (2*2)*(2-1)*ehp*ehp*(r**(-2*2)-1)
	             + (2*3)*(3+1)*thm*thm*(r**(2*3)-1) + (2*3)*(3-1)*thp*thp*(r**(-2*3)-1)
		     + (2*4)*(4+1)*shm*shm*(r**(2*4)-1) + (2*4)*(4-1)*shp*shp*(r**(-2*4)-1)
		     + (2*1)*(1+1)*dvm*dvm*(r**(2*1)-1) + (2*1)*(1-1)*dvp*dvp*(r**(-2*1)-1)
		     + (2*2)*(2+1)*evm*evm*(r**(2*2)-1) + (2*2)*(2-1)*evp*evp*(r**(-2*2)-1)
	             + (2*3)*(3+1)*tvm*tvm*(r**(2*3)-1) + (2*3)*(3-1)*tvp*tvp*(r**(-2*3)-1)
		     + (2*4)*(4+1)*svm*svm*(r**(2*4)-1) + (2*4)*(4-1)*svp*svp*(r**(-2*4)-1))/2$
%..............................................................................................

write "ff3v residual:";
ff3v - ff3v_target;

psi1  := psi1a + f1a * log(r)$

psi3p := ff3v/r$
psi3  := int(psi3p,r)$
psi3  := psi3a + psi3 - sub(r=1,psi3)$

%...............................................................................................
psi3_target := psi3a + f3a*log(r)
	       - f1a * ( - (1-dhpa+dhpa*dhpa + (2*2)*(2+1)*ehm*ehm + (2*2)*(2-1)*ehp*ehp
	       	       	   		     + (2*3)*(3+1)*thm*thm + (2*3)*(3-1)*thp*thp
					     + (2*4)*(4+1)*shm*shm + (2*4)*(4-1)*shp*shp
					     + (2*1)*(1+1)*dvm*dvm + (2*1)*(1-1)*dvp*dvp
					     + (2*2)*(2+1)*evm*evm + (2*2)*(2-1)*evp*evp
	       	       	   		     + (2*3)*(3+1)*tvm*tvm + (2*3)*(3-1)*tvp*tvp
					     + (2*4)*(4+1)*svm*svm + (2*4)*(4-1)*svp*svp)*log(r)
			 - dhpa*r*r*log(r) + r*r*log(r)*log(r)/2 + (1+dhpa*dhpa)*(r*r-1)/2
			 + (2+1)*ehm*ehm*(r**(2*2)-1) + (2-1)*ehp*ehp*(1-r**(-2*2))
			 + (3+1)*thm*thm*(r**(2*3)-1) + (3-1)*thp*thp*(1-r**(-2*3))
			 + (4+1)*shm*shm*(r**(2*4)-1) + (4-1)*shp*shp*(1-r**(-2*4))
			 + (1+1)*dvm*dvm*(r**(2*1)-1) + (1-1)*dvp*dvp*(1-r**(-2*1))
			 + (2+1)*evm*evm*(r**(2*2)-1) + (2-1)*evp*evp*(1-r**(-2*2))
			 + (3+1)*tvm*tvm*(r**(2*3)-1) + (3-1)*tvp*tvp*(1-r**(-2*3))
			 + (4+1)*svm*svm*(r**(2*4)-1) + (4-1)*svp*svp*(1-r**(-2*4)))/2$
% ..............................................................................................

write "psi3 residual:";
psi3 - psi3_target;

;bye;
