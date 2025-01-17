function [state,options,flag] = gaplotbestf1(options,state,flag)
drawnow;
if size(state.Score,2) > 1
    msg = getString(message('globaloptim:gaplotcommon:PlotFcnUnavailable','gaplotbestf'));
    title(msg,'interp','none');
    return;
end

app = get(findall(0,'Tag', 'MotorOptimisationAppTag'),'RunningAppInstance');
ga_ax = app.ga_UIAxes;
logarea = app.GAlogTextArea;
halt_b = app.ga_halt;

switch flag
    case 'init'
        hold(ga_ax,"on");
        set(ga_ax,'xlim',[0,options.MaxGenerations]);
        xlabel(ga_ax,'Generation','interp','none');
        ylabel(ga_ax,'Fitness value','interp','none');
        plotBest = plot(ga_ax,state.Generation,min(state.Score),'.k');
        set(plotBest,'Tag','gaplotbestf');
        plotMean = plot(ga_ax,state.Generation,meanf(state.Score),'.b');
        set(plotMean,'Tag','gaplotmean');
        title(ga_ax,['Best: ',' Mean: '],'interp','none');
        legend(ga_ax,'off');
        %         tp1 = sprintf('                              Best       Max        Stall');
        %         tp2 = sprintf('Generation  Func-count        f(x)     Constraint  Generations');
        tp = sprintf('\nSTARTING GA\n%s%s%s%s%s',pad('Generation',14,'both'), pad('Func-count',14,'both'), pad('f(x)',14,'both') ,pad('Constraint',14,'both'),pad('Stall Gens',14,'both'));
        lgval = get(logarea,'Value');
        lgval = cat(1,lgval, tp);
        set(logarea, 'Value', lgval);

    case 'iter'
        if(halt_b)
            state.StopFlag = 'STOP';
        end
        best = min(state.Score);
        m    = meanf(state.Score);
        plotBest = findobj(get(ga_ax,'Children'),'Tag','gaplotbestf');
        plotMean = findobj(get(ga_ax,'Children'),'Tag','gaplotmean');
        newX = [get(plotBest(1),'XData') state.Generation];
        newY = [get(plotBest(1),'YData') best];
        set(plotBest,'XData',newX, 'YData',newY);
        newY = [get(plotMean(1),'YData') m];
        set(plotMean,'XData',newX, 'YData',newY);
        set(get(ga_ax,'Title'),'String',sprintf('Best: %g Mean: %g',best,m));
        cineq = state.NonlinIneq;
        ceq = state.NonlinEq;
%         cnonlin = cat(1,ceq,cineq);
        maxConstr = 0;
        if ~isempty(cineq)
            maxConstr = max(cineq);
        end
        if length(ceq) > length(cineq)
            maxConstr = max([maxConstr;(ceq)]);
        end
        gen = string(state.Generation);
        feval = string(state.FunEval);
        bst = string(state.Best(end));
        mxcstr = string(abs(maxConstr));
        stllgen = string(state.Generation - state.LastImprovement);
        tp = sprintf('%s%s%s%s%s', pad(gen,14,'both'), pad(feval,14,'both'), pad(bst,14,'both') ,pad(mxcstr,14,'both'),pad(stllgen,14,'both'));
        lgval = get(logarea,'Value');
        lgval = cat(1,lgval, tp);
        set(logarea, 'Value', lgval);
%         pause(0.15);
        drawnow;
    case 'done'
        LegnD = legend(ga_ax,'Best fitness','Mean fitness');
        set(LegnD,'FontSize',8);
        hold(ga_ax,"off");
end
end
%------------------------------------------------
function m = meanf(x)
nans = isnan(x);
x(nans) = 0;
n = sum(~nans);
n(n==0) = NaN; % prevent divideByZero warnings
% Sum up non-NaNs, and divide by the number of non-NaNs.
m = sum(x) ./ n;
end



