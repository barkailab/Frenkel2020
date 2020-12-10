function [rt,p]=calculate_RT_sigmafit(data,meta,norm_vec,tp,varargin)
ip=inputParser;
ip.addParameter('th',0.5);
ip.addParameter('fit','sigma');
ip.addParameter('ws',2);
ip.addParameter('gb',true(size(meta,1),1));
ip.addParameter('ylog',true);
ip.addParameter('disTh',10);

ip.parse(varargin{:});

th=ip.Results.th;
fit = find(contains({'sigma','linear','robust'},ip.Results.fit));% sigma =1; linear =2; robust =3
%final=2.^norm_vec(end)-1
%sigfunc=@(A,x)log2(1+(final./(1+exp(-(x-A(2))/+A(1)))));
sigfunc=@(A,x)norm_vec(end)./(1+exp(-(x-A(2))/+A(1)));
weights=(filter2([1 1 1],norm_vec>0 & norm_vec<0.95*norm_vec(end))>0)*3+1;
p=nan(size(data,1),2);
rt=nan(size(data,1),numel(th));
s=nan(size(data,1),1);
e=nan(size(data,1),1);
gb_fit=conv2(ones(3,1),1,ip.Results.gb,'same')>0;
for i=find(gb_fit)'
    sel_i=i-ip.Results.ws:i+ip.Results.ws;
    sel_i=sel_i(sel_i>0 & sel_i<size(data,1));
    if fit ==1
        x_i=acol(repmat(tp,numel(sel_i),1));
        if ip.Results.ylog
            y_i=acol(data(sel_i,:)+norm_vec);
        else
            y_i=2.^acol(data(sel_i,:)+norm_vec)-1;
        end
        w_i=acol(repmat(weights,numel(sel_i),1));
        oki=~isnan(y_i);
        if sum(oki)>2*numel(tp) & sum(x_i(oki)>.5)>numel(sel_i)
            score_func=@(A) sum(((sigfunc(A,x_i(oki))-y_i(oki)).^2).*w_i(oki))/sum(w_i(oki));
            [p(i,:),s(i),e(i)]= fminsearch(score_func,[5 25],optimset('MaxFunEvals',5000));
            for c=1:numel(th)
                rt(i,c)=min([find(sigfunc(p(i,:),0.01:.01:50)>th(c),1)*.01,NaN]);
            end
        end
    elseif fit == 2 | fit ==3
        sel_i=sel_i(all(~isnan(data(sel_i,:)),2));        
        tp_i=find(movmean(mean(data(sel_i,:)+norm_vec,1)>=mean(th),5)>0.5,1); %select best timepont for linear fit
        tp_x(i)=max([NaN,tp_i]);
        tp_i=tp_i-2:tp_i+2;
        tp_i=tp_i(tp_i>0 & tp_i<numel(tp));
        if numel(tp_i)>=3 
            fit_i=nan(1,numel(tp_i));
            if range(meta(sel_i,2))>0 % calculate dna content for fit either median or based 
                for t=1:numel(tp_i)
                    p_t=polyfit(meta(sel_i,2),data(sel_i,tp_i(t)),1);
                    fit_i(t)=meta(i,2)*p_t(1)+p_t(2);
                end
            else
                fit_i=median(data(sel_i,tp_i),1);
            end
            if ip.Results.ylog
                fit_i=fit_i+norm_vec(tp_i);
            else
                fit_i=2.^(fit_i+norm_vec(tp_i))-1;
            end
            if fit == 3
                pi=robustfit(tp(tp_i),fit_i);
                p(i,:)=[pi(2) pi(1)];
            else
                p(i,:)=polyfit(tp(tp_i),fit_i,1); 
            end
            for c=1:numel(th)
                rt(i,c)=(th(c)-p(i,2))/p(i,1);
                if min(abs(rt(i,c)-tp(tp_i)))>ip.Results.disTh | p(i,1)<0
                    rt(i,c)=NaN;
                end
            end
        end
    end
end
rt(rt>=tp(end))=NaN;
end