function [pu_m,pu_c,a] = piecewiseUnder(x,cdf,usererror,pieces)
% This function provides the slope and y-intercepts for the given points.


    maxerror = [];
    N = length(x);
    xpoints(:,1) = x;
    ypoints(:,1) = cdf; 
    continueCondition = 1;
    pu_m = [];
    pu_c = [];
    
    while continueCondition == 1
        
        xeval = triu(ones(length(xpoints)-1,length(xpoints))).*xpoints';
        yeval = triu(ones(length(ypoints)-1,length(ypoints))).*ypoints';
        
        xdiff = xpoints(end)-xpoints(1:end-1);
        ydiff = ypoints(end)-ypoints(1:end-1);
        m(1,:) = ydiff./xdiff;
        c(1,:) = ypoints(end) - m'.*xpoints(end);

        y = m.*xeval'+c;
        y = tril(y);
        
        error = yeval'- y;

        % Remove those slopes and intercepts which have negative error: 
        
        poserror = all(error >= 0,1);
        maxerror = max(abs(error))<=usererror;
        totmet = poserror+maxerror;
        
        metric = find(totmet==2);
%         hold on, plot(xeval(metric(1),:),yeval(metric(1),:)'-y(:,metric(1)))
        
        if length(metric) == 1 || isempty(metric) == 1 ||  length(pu_m) >= pieces
            break;
        else
            xpoints = xpoints(1:metric(1));
            ypoints = ypoints(1:metric(1));
            pu_m = [pu_m m(metric(1))];
            pu_c = [pu_c c(metric(1))];
            a(1) = xpoints(metric(1));
            a(2) = ypoints(metric(1));
            xeval = [];
            yeval = [];
            xdiff = [];
            ydiff = [];
            m = [];
            c = [];
            y = [];
            error = [];
            poserror = [];
            tolerror = [];
            totmet = [];
            metric = [];
        end
        
        
    end

end

