function law=makeHDivBHMaterial(table)
% Bounded monotone soft-iron inverse BH, CPU production PCHIP contract.
validateattributes(table,{'double'},{'2d','real','finite'});
assert(size(table,2)==2 && size(table,1)>=3 && size(table,1)<=4096 && all(table(1,:)==0));
h=table(:,1);b=table(:,2);mu0=4e-7*pi;mag=b/mu0-h;
assert(all(diff(h)>0) && all(diff(b)>0) && all(diff(mag)>0),'Strictly increasing origin-anchored magnetization required');
forward=pchip(h,mag);cb=mu0*forward.coefs;cb(:,3)=cb(:,3)+mu0;cb(:,4)=cb(:,4)+mu0*h(1:end-1);forward=mkpp(h,cb);
hs=[0 logspace(-2,log10(max(h(end),1)),800)];ms=ppval(forward,hs)/mu0-hs;
ms(hs>h(end))=mag(end);keep=false(size(ms));last=-inf;
for k=1:numel(ms),if ms(k)>last+1e-12*max(1,mag(end)),keep(k)=true;last=ms(k);end,end
assert(nnz(keep)>=2,'Inverse BH grid is not monotone');
pp=pchip(ms(keep),hs(keep));knots=pp.breaks(:);c=pp.coefs;dx=diff(knots);
areas=c(:,1).*dx.^4/4+c(:,2).*dx.^3/3+c(:,3).*dx.^2/2+c(:,4).*dx;
base=[0;cumsum(areas)];cap=mag(end)*(1-1e-9);
assert(cap>=0 && cap<=knots(end));
idx=find(knots<=cap,1,'last');idx=min(idx,size(c,1));x=cap-knots(idx);
hlast=polyval(c(idx,:),x);nu0=c(1,3);
assert(isfinite(nu0) && nu0>0 && all(isfinite(c(:))));
law=struct('schema_version',1,'knots',knots,'coefficients',c,'integral_base',base, ...
    'cap',cap,'barrier_slope',1000,'last_field',hlast,'initial_reluctivity',nu0, ...
    'table',table,'units','H and M A/m; B tesla; coenergy A^2/m^2 (common mu0 factor omitted)');
end
