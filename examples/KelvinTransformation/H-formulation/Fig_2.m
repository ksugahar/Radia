clear all;
close all;

set(gcf,'Units','Pixels','Position',[100,30,1200,1170],'PaperSize',[1200/96*2.54,1170/96*2.54]);
for m = 1:2
for n = 1:2
	a(m,n) = axes('Units','Pixels','Position',[155+(n-1)*600,135+(2-m)*600,425,425],'FontName','Times New Roman','FontSize',36);
	hold on;	box on;	grid on; grid minor;
	axis equal;
	xlabel('{\it x} (m)');
	ylabel('{\it z} (m)');
	set(gca,'XLim',[-1.1, 1.1], 'XTick', [-1.0:1.00:1.0]);
	set(gca,'YLim',[-1.1, 1.1]);
end
end

theta = linspace(0,2*pi,101);

for m = 1:2
for n = 1:2
	axes(a(m,n))
		if isequal([m,n], [1,1])
			load('3D_quadrupole_with_Kelvin.mat');
			Hx = Hx_analytical;
			Hz = Hz_analytical;
			H = sqrt(Hx_analytical.^2+Hz_analytical.^2);
			text(0.75,-1.6,'(a)','FontName','Times New Roman','FontSize', 36);
			h1(1) = plot(0.5*cos(theta), 0.5*sin(theta), 'Color', 'k', 'LineWidth', 0.5);
			h1(2) = text(0,0,'{\it\mu}=100','HorizontalAlignment','center', 'FontName','Times New Roman','FontSize', 36);
		elseif isequal([m,n], [1,2])
			load('3D_quadrupole.mat');
			H = sqrt(Hx.^2+Hz.^2);
			text(0.75,-1.6,'(b)','FontName','Times New Roman','FontSize', 36);
			h1(1) = plot(0.5*cos(theta), 0.5*sin(theta), 'Color', 'k', 'LineWidth', 0.5);
			h1(2) = text(0,0,'{\it\mu}=100','HorizontalAlignment','center', 'FontName','Times New Roman','FontSize', 36);
		elseif isequal([m,n], [2,1])
			load('3D_quadrupole_with_Kelvin.mat');
			H = sqrt(Hx.^2+Hz.^2);
			text(0.75,-1.6,'(c)','FontName','Times New Roman','FontSize', 36);
			h1(1) = plot(0.5*cos(theta), 0.5*sin(theta), 'Color', 'k', 'LineWidth', 0.5);
			h1(2) = text(0,0,'{\it\mu}=100','HorizontalAlignment','center', 'FontName','Times New Roman','FontSize', 36);
		elseif isequal([m,n], [2,2])
			load('3D_quadrupole_with_Kelvin.mat');
			xx = xx_ext-3;
			zz = zz_ext;
			Hx = Hx_ext;
			Hz = Hz_ext;
			H = sqrt(Hx.^2+Hz.^2);
			text(0.75,-1.6,'(d)','FontName','Times New Roman','FontSize', 36);
		end
		[c,h] = contourf(xx, zz, H, linspace(0, 1.8, 30), 'LineColor', 'none');
		h = streamslice(xx, zz, Hx, Hz);
	set(h, 'Color', 'w', 'LineWidth', 0.1);

		uistack(h1, 'top');

		set(h, 'Color', 'w', 'LineWidth', 0.1);
end
end

exportgraphics(gcf, sprintf('%s.pdf', mfilename));
