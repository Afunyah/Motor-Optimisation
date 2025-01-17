function stop = optimplotfval1(x,optimValues,state)
drawnow;
app = get(findall(0,'Tag', 'MotorOptimisationAppTag'),'RunningAppInstance');
logarea = app.GAlogTextArea;
halt_b = app.fmin_halt;
fmin_ax = app.fmin_UIAxes;

stop = false;
switch state
    case 'init'
        tp = sprintf('\nSWITCHING TO FMINCON\n%s%s%s%s%s',pad('Iter',14,'both'), pad('Func-count',14,'both'), pad('f(x)',14,'both'), pad('Constraint',14,'both'),pad('Optimality',14,'both'));
        lgval = get(logarea,'Value');
        lgval = cat(1,lgval, tp);
        set(logarea, 'Value', lgval);
    case 'iter'
        if(halt_b)
            stop = true;
        end
        if isfield(optimValues,'fval')
            if isscalar(optimValues.fval)
                plotscalar(optimValues.iteration,optimValues.fval,fmin_ax);
                itr = string(optimValues.iteration);
                fcnt = string(optimValues.funccount);
                feval = string(optimValues.fval);
                cnstr = string(optimValues.constrviolation);
                optm = string(optimValues.firstorderopt);
                tp = sprintf('%s%s%s%s%s', pad(itr,14,'both'), pad(fcnt,14,'both'), pad(feval,14,'both'), pad(cnstr,14,'both'),pad(optm,14,'both'));
                lgval = get(logarea,'Value');
                lgval = cat(1,lgval, tp);
                set(logarea, 'Value', lgval);
            end
        end
    case 'done'
end

function plotscalar(iteration,fval,fmin_ax)



if iteration == 0
    plotfval = plot(fmin_ax,iteration,fval,'kd','MarkerFaceColor',[1 0 1]);
    title(fmin_ax,getString(message('MATLAB:optimfun:funfun:optimplots:TitleCurrentFunctionValue',sprintf('%g',fval))),'interp','none');
    xlabel(fmin_ax,getString(message('MATLAB:optimfun:funfun:optimplots:LabelIteration')),'interp','none');
    set(plotfval,'Tag','optimplotfval');
    ylabel(fmin_ax,getString(message('MATLAB:optimfun:funfun:optimplots:LabelFunctionValue')),'interp','none')
    drawnow;
else
    plotfval = findobj(get(fmin_ax,'Children'),'Tag','optimplotfval');
    newX = [get(plotfval,'Xdata') iteration];
    newY = [get(plotfval,'Ydata') fval];
    set(plotfval,'Xdata',newX, 'Ydata',newY);
    set(get(fmin_ax,'Title'),'String',getString(message('MATLAB:optimfun:funfun:optimplots:TitleCurrentFunctionValue',sprintf('%g',fval))));
    drawnow;
end
