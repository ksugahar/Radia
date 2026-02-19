clear all;
close all;

set(gcf,'Units','Pixels','Position',[100,30,1200,500],'PaperSize',[1200/96*2.54,500/96*2.54]);
for n = 1:2
	a(n) = axes('Units','Pixels','Position',[155+(n-1)*600,135,425,345],'FontName','Times New Roman','FontSize',36);
	hold on;	box on;	grid on; grid minor;
	xlabel('{\it x} (m)');
	ylabel('{\it H_z} (A/m)');
	set(gca,'XLim',[-1.1, 1.1], 'XTick', [-1.0:1.00:1.0]);
end

theta = linspace(0,2*pi,101);

for n = 1:2
	axes(a(n))
		switch n
		case 1
			load('3D_dipole_with_Kelvin.mat');
			x = xx(111,:);
			Hz = Hz_analytical(111,:);
			Hz(find(abs(x)>=0.5)) = - Hz(find(abs(x)>=0.5));
			plot(x, Hz, 'k-', 'LineWidth', 1.0);

			load('3D_dipole.mat');
			x = xx(111,:);
			Hz = Hz(111,:);
			plot(x, Hz, 'r-');

			load('3D_dipole_with_Kelvin.mat');
			x = xx(111,:);
			Hz = Hz(111,:);
			plot(x, Hz, 'c--', 'LineWidth', 1.5);

			set(gca,'YLim',[-1.1, 0.0]);
			text(0.75,-1.4,'(a)','FontName','Times New Roman','FontSize', 36);
			h = legend({'Analytical', 'Finite domain' 'Kelvin trans.'},'location','best','box','off', 'FontName','Times New Roman','FontSize', 28);
		case 2
			load('3D_quadrupole_with_Kelvin.mat');
			x = xx(111,:);
			Hz = Hz_analytical(111,:);
			Hz(find(abs(x)>=0.5)) = Hz(find(abs(x)>=0.5));
			plot(x, Hz, 'k-', 'LineWidth', 1.0);

			load('3D_quadrupole.mat');
			x = xx(111,:);
			Hz = Hz(111,:);
			plot(x, Hz, 'r-');

			load('3D_quadrupole_with_Kelvin.mat');
			x = xx(111,:);
			Hz = Hz(111,:);
			plot(x, Hz, 'c--', 'LineWidth', 1.5);

			set(gca,'YLim',[-0.60, 0.60]);
			text(0.75,-0.93,'(b)','FontName','Times New Roman','FontSize', 36);
			h = legend({'Analytical', 'Finite domain' 'Kelvin trans.'},'location','southeast','box','off', 'FontName','Times New Roman','FontSize', 28);
		end
end

exportgraphics(gcf, sprintf('%s.pdf', mfilename));
