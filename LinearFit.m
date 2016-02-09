function [ tconsts, Anoughts ] = LinearFit( matFile, Tmin, Tmax )
%load('atriK1.mat')
%load(matFile);

t = load(matFile,'T');  % Function output form of LOAD
t = t.('T');
% index from 251 to 1000, for 15 < T < 60
% fit only for T > 15
Tindex = (t > Tmin) & (t < Tmax);
xdata = t(Tindex);

configs = {'o1','o2','o3','o4','o5','oA'};
tconsts = zeros(1,size(configs,2));
Anoughts = tconsts;

%rawdata = zeros(size(t,2),size(configs,2));
for config = 1:size(configs,2)
	con = load(matFile,configs{config});  % Function output form of LOAD
	con = con.(configs{config});
	%rawdata(:,config) = con;

	% calculate natural log, so fit straight line
	ydata = log(con(Tindex));

	% fitting - http://www.mathworks.com/help/matlab/data_analysis/programmatic-fitting.html#f1-7407
	[decay,ErrorEst] = polyfit(xdata,ydata,1);
	tconsts(config) = decay(1);
	Anoughts(config) = decay(2);
	decay_fit = polyval(decay,xdata,ErrorEst);

	% plot data and fit
	h = figure;
	plot(xdata,decay_fit,'-',xdata,ydata,'+');
	title([configs{config},'-t1-',num2str(Tmin),'-t2-',num2str(Tmax),'-grad=',num2str(decay(1))])
	legend('Linear Model','Data','Location','NorthEast');
	axis([Tmin,Tmax,-Inf,Inf])
	xlabel('Integration time');
	ylabel('log correlation');
	print(h,'-dpng',[matFile,'-',configs{config},'-t1-',num2str(Tmin),'-t2-',num2str(Tmax),'-LOGplot.png']);

	% residuals
	res = ydata - decay_fit;
	%figure, plot(xdata,res,'+')
	%title('Residuals for the Linear Model')
end

end
