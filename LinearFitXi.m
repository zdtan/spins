function [ xis, intercepts ] = LinearFitXi( matFile, cutoff )
%load('apicolK1.mat')
%load(matFile);

x = load(matFile,'sites');  % Function output form of LOAD
x = x.('sites');

%% note that x starts from 1. we want to start from 0.
% want to collapse by averaging over periodic boundary conditions
Nsites = size(x,2);
x2 = [0:floor(Nsites/2)];	% starts from 0

configs = {'c1','c2','c3','c4','c5','cA'};
xis = zeros(1,size(configs,2));
intercepts = xis;
for config = 1:size(configs,2)
	con = load(matFile,configs{config});  % Function output form of LOAD
	con = con.(configs{config});
	con2 = zeros(1,size(x2,2));
	for pos = 1+x2	% starts from 1
		con2(pos) = con(pos);
		if pos > 1
			con2(pos) = ( con2(pos) + con(end+2-pos) )/2;
		end
	end

	% cutoff for data
	sIndex = con2 > cutoff;
	xdata = x2(sIndex);

	% calculate natural log, so fit straight line
	ydata = log(con2(sIndex));

	% fitting - http://www.mathworks.com/help/matlab/data_analysis/programmatic-fitting.html#f1-7407
	[decay,ErrorEst] = polyfit(xdata,ydata,1);
	xis(config) = decay(1);
	intercepts(config) = decay(2);
	decay_fit = polyval(decay,xdata,ErrorEst);

	% plot data and fit
	h = figure;
	plot(xdata,decay_fit,'-',xdata,ydata,'+');
	title([configs{config},'-cut',num2str(cutoff),'-grad=',num2str(decay(1))])
	legend('Linear Model','Data','Location','NorthEast');
	axis([0,floor(Nsites/2),-Inf,0])
	xlabel('spin separation along apical column');
	ylabel('log correlation');
	print(h,'-dpng',[matFile,'-',configs{config},'-cut',num2str(cutoff),'-LOGplot.png']);

	% residuals
	res = ydata - decay_fit;
	%figure, plot(xdata,res,'+')
	%title('Residuals for the Linear Model')
end

end
