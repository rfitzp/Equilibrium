% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to perform matching between near vacuum and far vacuum solutions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linelength 120$

out "Matching.out"$

write "<< Output from Matching.red >>";

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Poloidal flux from near-vacuum region
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi0 := psi1a + f1a * log(r)$

psi2 := psi3a + f3a*log(r)
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
			 
psi := psi0 + eps**2*psi2$

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Poloidal flux from far-vacuum region
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi0 := 1$

xi1 := r*cos(eta)/2 + cos(pi)*dh*cos(eta)/r + cos(2*pi)*eh*cos(2*eta)/r + cos(3*pi)*th*cos(3*eta)/r + cos(4*pi)*sh*cos(4*eta)/r
	            - cos(pi)*dv*sin(eta)/r - cos(2*pi)*ev*sin(2*eta)/r - cos(3*pi)*tv*sin(3*eta)/r - cos(4*pi)*sv*sin(4*eta)/r$
		    
xi2 := r*r/4 + dh*dh/2/r/r + eh*eh/2/r/r + th*th/2/r/r + sh*sh/2/r/r
	     + dv*dv/2/r/r + ev*ev/2/r/r + tv*tv/2/r/r + sv*sv/2/r/r$

xi := xi0 + eps*xi1 + eps**2*xi2$

xir := xi/r$
rxi := r * av2( inv(xi,eps), eps, eta)$

fac     := 1 - eps*rxi*cos(eta)$
factor1 := av2( pow2(fac,eps,-1,2), eps, eta)$

lnxi0 := 0$

lnxi1 := r*cos(eta)/2 + cos(pi)*dh*cos(eta)/r + cos(2*pi)*eh*cos(2*eta)/r + cos(3*pi)*th*cos(3*eta)/r + cos(4*pi)*sh*cos(4*eta)/r
                      - cos(pi)*dv*sin(eta)/r - cos(2*pi)*ev*sin(2*eta)/r - cos(3*pi)*tv*sin(3*eta)/r - cos(4*pi)*sv*sin(4*eta)/r$

lnxi2 := r*r/4 + dh*dh/2/r/r + eh*eh/2/r/r + th*th/2/r/r + sh*sh/2/r/r
	     	      + dv*dv/2/r/r + ev*ev/2/r/r + tv*tv/2/r/r + sv*sv/2/r/r$

lnxi := lnxi0 + eps*lnxi1 + eps**2*lnxi2$

factor2 :=   pc0 * (1 - l8ea/2 + log(r)/2 - lnxi/2) - pc0 * (eps**2)*(rxi**2) * (1 - l8ea + log(r) - lnxi)/32
	   + qc0 * (1 - (eps***2)*(rxi**2)/16)
	   + pc1 * eps*xir     *cos(eta)   + qc1 * eps*rxi     *cos(eta)
	   + pc2 * eps*(xir**2)*cos(2*eta) + qc2 * eps*(rxi**2)*cos(2*eta)
	   + pc3 * eps*(xir**3)*cos(3*eta) + qc3 * eps*(rxi**3)*cos(3*eta)
	   + pc4 * eps*(xir**4)*cos(4*eta) + qc4 * eps*(rxi**4)*cos(4*eta)
	   + ps1 * eps*xir     *sin(eta)   + qs1 * eps*rxi     *sin(eta)
	   + ps2 * eps*(xir**2)*sin(2*eta) + qs2 * eps*(rxi**2)*sin(2*eta)
	   + ps3 * eps*(xir**3)*sin(3*eta) + qs3 * eps*(rxi**3)*sin(3*eta)
	   + ps4 * eps*(xir**4)*sin(4*eta) + qs4 * eps*(rxi**4)*sin(4*eta)$

factor2 = av2(factor2,eps,eta)$

psivac := av2(factor1*factor2,eps,eta)$

psivac0   := coeffn(psivac,eps,0)$
psivac1   := coeffn(psivac,eps,1)$
psivac1c0 := fcos(psivac1,eta,0)$
psivac1c1 := fcos(psivac1,eta,1)$
psivac1c2 := fcos(psivac1,eta,2)$
psivac1c3 := fcos(psivac1,eta,3)$
psivac1c4 := fcos(psivac1,eta,4)$
psivac1s1 := fsin(psivac1,eta,1)$
psivac1s2 := fsin(psivac1,eta,2)$
psivac1s3 := fsin(psivac1,eta,3)$
psivac1s4 := fsin(psivac1,eta,4)$
psivac2   := coeffn(psivac,eps,2)$

%...................................................................................
psivac0_target   := pc0 * (1 - l8ea/2 + log(r)/2) + qc0$

psivac1c1_target := pc0 * (r - r*l8ea + r*log(r))/4 + qc0*r/2
                    - pc0*cos(pi)  *dh/2/r + pc1/r       + qc1*r$
psivac1c2_target := - pc0*cos(2*pi)*eh/2/r + pc2/r/r     + qc2*r*r$
psivac1c3_target := - pc0*cos(3*pi)*th/2/r + pc3/r/r/r   + qc3*r*r*r$
psivac1c4_target := - pc0*cos(4*pi)*sh/2/r + pc4/r/r/r/r + qc4*r*r*r*r$

psivac1s1_target :=   pc0*cos(pi)  *dv/2/r + ps1/r       + qs1*r$
psivac1s2_target :=   pc0*cos(2*pi)*ev/2/r + ps2/r/r     + qs2*r*r$
psivac1s3_target :=   pc0*cos(3*pi)*tv/2/r + ps3/r/r/r   + qs3*r*r*r$
psivac1s4_target :=   pc0*cos(4*pi)*sv/2/r + ps4/r/r/r/r + qs4*r*r*r*r$

psivac2_target   := pc0 * (-5*r*r/32 - (-3/8 + l8ea/8) * dh + dh*log(r)/8
                          - dh*dh/4/r/r - eh*eh/4/r/r - th*th/4/r/r - sh*sh/4/r/r
			  - dv*dv/4/r/r - ev*ev/4/r/r - tv*tv/4/r/r - sv*sv/4/r/r)
		   + qc0 * dh/4
		   + pc1 * 1/2 
		   + pc1*cos(pi)  *1*dh/r      /2/r - qc1*cos(pi)  *1*dh*r      /2/r
		   + pc2*cos(2*pi)*2*eh/r/r    /2/r - qc2*cos(2*pi)*2*eh*r*r    /2/r
		   + pc3*cos(3*pi)*3*th/r/r/r  /2/r - qc3*cos(3*pi)*3*th*r*r*r  /2/r
		   + pc4*cos(4*pi)*4*sh/r/r/r/r/2/r - qc4*cos(4*pi)*4*sh*r*r*r*r/2/r
		   - ps1*cos(pi)  *1*dv/r      /2/r + qs1*cos(pi)  *1*dv*r      /2/r
		   - ps2*cos(2*pi)*2*ev/r/r    /2/r + qs2*cos(2*pi)*2*ev*r*r    /2/r
		   - ps3*cos(3*pi)*3*tv/r/r/r  /2/r + qs3*cos(3*pi)*3*tv*r*r*r  /2/r
		   - ps4*cos(4*pi)*4*sv/r/r/r/r/2/r + qs4*cos(4*pi)*4*sv*r*r*r*r/2/r$
%....................................................................................
		  
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

% %%%%%%%%%%%%%%%%%%%%%
% Zeroth order matching
% %%%%%%%%%%%%%%%%%%%%%

dh := dha - r*r*log(r)/2 + (2*dhpa + 1)*(r*r - 1)/4$
eh := ehm * r**(2+1) - ehp * r**(1-2)$
th := thm * r**(3+1) - thp * r**(1-3)$
sh := shm * r**(4+1) - shp * r**(1-4)$
dv := dvm * r**(1+1) - dvp * r**(1-1)$
ev := evm * r**(2+1) - evp * r**(1-2)$
tv := tvm * r**(3+1) - tvp * r**(1-3)$
sv := svm * r**(4+1) - svp * r**(1-4)$

%................................
pc0 := 2*f1a$
qc0 := psi1a - 2*f1a*(1 - l8ea/2)$
%................................

residulal0 := psi0 - psivac0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order cos(eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%............................................
pc1 := f1a * (1/4 - dha + dhpa/2)$
qc1 := - qc0/2 + f1a * (l8ea - 3/2 - dhpa)/2$
%............................................

residual1c1 := psivac1c1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order cos(2*eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..........................
pc2 := - f1a*cos(2*pi)*ehp$
qc2 :=   f1a*cos(2*pi)*ehm$ 
%..........................

residual1c2 := psivac1c2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order cos(3*eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..........................
pc3 := - f1a*cos(3*pi)*thp$
qc3 :=   f1a*cos(3*pi)*thm$ 
%..........................

residual1c3 := psivac1c3;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order cos(4*eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..........................
pc4 := - f1a*cos(4*pi)*shp$
qc4 :=   f1a*cos(4*pi)*shm$ 
%..........................

residual1c4 := psivac1c4;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order sin(eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%........................
ps1 :=   f1a*cos(pi)*dvp$
qs1 := - f1a*cos(pi)*dvm$ 
%........................

residual1s1 := psivac1s1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order sin(2*eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..........................
ps2 :=   f1a*cos(2*pi)*evp$
qs2 := - f1a*cos(2*pi)*evm$ 
%..........................

residual1s2 := psivac1s2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order sin(3*eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..........................
ps3 :=   f1a*cos(3*pi)*tvp$
qs3 := - f1a*cos(3*pi)*tvm$ 
%..........................

residual1s3 := psivac1s3;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order sin(4*eta) matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..........................
ps4 :=   f1a*cos(4*pi)*svp$
qs4 := - f1a*cos(4*pi)*svm$ 
%..........................

residual1s4 := psivac1s4;

% %%%%%%%%%%%%%%%%%%%%%
% Second order matching
% %%%%%%%%%%%%%%%%%%%%%

residual2 := psi2 - psivac2$

residual20 := coeffn(residual2, log(r), 0)$
residual21 := coeffn(residual2, log(r), 1)$

residual2 := residual2 - residual20 - residual21*log(r);

%..................................................................................
residual20_target := psi3a + f1a * (3/16 + dha/4 - dhpa/4 + dha*dhpa/2
		     	     	    - (2-1)*ehp*ehp/2 - ehp*ehm + (2+1)*ehm*ehm/2
				    - (3-1)*thp*thp/2 - thp*thm + (3+1)*thm*thm/2
				    - (4-1)*shp*shp/2 - shp*shm + (4+1)*shm*shm/2
				    - (1-1)*dvp*dvp/2 - dvp*dvm + (1+1)*dvm*dvm/2
				    - (2-1)*evp*evp/2 - evp*evm + (2+1)*evm*evm/2
				    - (3-1)*tvp*tvp/2 - tvp*tvm + (3+1)*tvm*tvm/2
				    - (4-1)*svp*svp/2 - svp*svm + (4+1)*svm*svm/2)$
residual21_target := f3a   + f1a * (5/8 - dha/2 - dhpa/4 + dhpa*dhpa/2
		     	     	    + 2*(2-1)*ehp*ehp + 2*(2+1)*ehm*ehm
				    + 3*(3-1)*thp*thp + 3*(3+1)*thm*thm
				    + 4*(4-1)*shp*shp + 4*(4+1)*shm*shm
				    + 1*(1-1)*dvp*dvp + 1*(1+1)*dvm*dvm
				    + 2*(2-1)*evp*evp + 2*(2+1)*evm*evm
				    + 3*(3-1)*tvp*tvp + 3*(3+1)*tvm*tvm
				    + 4*(4-1)*svp*svp + 4*(4+1)*svm*svm)$
%..................................................................................				    

residual20 := residual20 - residual20_target;
residual21 := residual21 - residual21_target;

out t$
;bye;