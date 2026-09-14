clear all;
close all;

options = optimset('TolFun',1e-15,'TolX',1e-15,'MaxIter',10000,'MaxFunEvals',20000,'Display','none');

set(gcf,'Units','Pixels','Position',[2560,30,1330,590]);
a1 = axes('Units','Pixels','Position',[100,80,500,500],'FontName','Times New Roman','FontSize',20);
	hold on;	box on;	grid on;	grid minor;

	S = 0.9;
	Ms = 1.91;
	Hc = 500;
	a = 1.0;
	HMax = linspace(0,4000,41);
	HMax = HMax(2:end);
	dH = HMax(1);
	BH = [];
	for Hmax = HMax
		CF = @(a) match_Hmax(Hmax, a, S, Ms, Hc);
		[a, fval] = fminsearch(CF,[a],options);
		disp(sprintf('error = %.2e', fval));
		H = -Hmax:dH:Hmax;
		M = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
		a = -a;
		M = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
		a = -a;
		BH(length(BH)+1).H = H(end:-1:1);
		BH(length(BH)+0).B = M(end:-1:1);
	end

	BH = BH(end:-1:1);
	for n = 1:length(BH)
		h(1) = plot( BH(n).H, BH(n).B, 'r-');
		h(2) = plot(-BH(n).H,-BH(n).B, 'b-');
	end
	save(sprintf('H_input.mat'),'BH','dH','HMax');

	legend(h, {sprintf('Descending'), sprintf('Ascending\n ({\\Delta}{\\itH}=%.1f A/m)',dH)}, 'box','off', 'location','northwest');
	set(gca,'XLim',[-2000, 2000]);
	xlabel('{\itH} (A/m)');
	ylabel('{\itB} (T)');

a2 = axes('Units','Pixels','Position',[800,80,500,500],'FontName','Times New Roman','FontSize',20);
	hold on;	box on;	grid on;	grid minor;

	BMax = linspace(0,1.9,21);
	BMax = BMax(2:end);
	dB = BMax(1);
	BH = [];
	for Bmax = BMax
		Hmax = 100;
		CF = @(x) match_Mmax(Bmax, x(1), x(2), S, Ms, Hc);
		[x, fval] = fminsearch(CF,[Hmax; a], options);
		disp(sprintf('error = %.2e', fval));
		Hmax = x(1);
		a = x(2);
		H = linspace(-Hmax, Hmax, 2000);
		M = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
		a = -a;
		M = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
		a = -a;
		B = Bmax:-dB:-Bmax;
		H = interp1(M,H,B,'makima');
		BH(length(BH)+1).H = H;
		BH(length(BH)+0).B = B;
	end

	BH = BH(end:-1:1);
	for n = 1:length(BH)
		h(1) = plot( BH(n).H, BH(n).B, 'r-');
		h(2) = plot(-BH(n).H,-BH(n).B, 'b-');
	end
	save(sprintf('B_input.mat'),'BH','dB','BMax');

	legend(h, {sprintf('Descending'), sprintf('Ascending\n ({\\Delta}{\\itB}=%.3f T) ',dB)}, 'box','off', 'location','northwest');
	set(gca,'XLim',[-2000, 2000]);
	xlabel('{\itH} (A/m)');
	ylabel('{\itM} (T)');

print('-dpng',sprintf('%s.png',mfilename));

function dM = match_Hmax(H, a, S, Ms, Hc);
	M1 = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
	a = -a;
	M2 = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
	dM = (M1-M2)^2;
end

function dM = match_Mmax(M, H, a, S, Ms, Hc);
	M1 = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
	a = -a;
	M2 = Ms*(sign(a) - a*(1+tanh((Hc-sign(a)*H)/Hc*atanh(S))));
	dM = (M - M1)^2 + (M - M2)^2;;
end
