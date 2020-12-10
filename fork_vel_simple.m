function out=fork_vel_simple(rt,meta,varargin)
ip=inputParser;
ip.addParameter('axis',[]);
ip.addParameter('ulim',11000);
ip.addParameter('llim',3000);
ip.addParameter('adjOri',false);
ip.addParameter('rtth',22);
ip.addParameter('fire',[]);
ip.addParameter('stop',false);
ip.addParameter('filterOri',false);
ip.addParameter('oldOris',false);

ip.parse(varargin{:})
ulim=ip.Results.ulim;
llim=ip.Results.llim;
sels=~isnan(meta(:,4));
if ip.Results.oldOris
    load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/DNAreplication/pubdata/oris.mat', 'oridb')
    GP=load('group_imp.mat');
    chr_idx=[0 cumsum(GP.chr_len)];
    oridb(:,8)=chr_idx(oridb(:,1))'+mean(oridb(:,2:3),2);
    [~,idx]=min(abs(GP.oris.loc-oridb(:,8)'));
    oris=GP.oris.loc(unique(idx));
    clear dis idx oridb GP chr_idx
    dismat=meta(:,1)-oris';
    [~,idx]=min(abs(dismat),[],2);
    meta(:,3:4)=[dismat(sub2ind(size(dismat),1:24315,idx'))' idx];
    clear dismat idx oris
end
min_dis=accumarray(meta(sels,4),abs(meta(sels,3)),[],@min);
oris=find(abs(meta(sels,3))==min_dis(meta(sels,4)));

clear sels min_dis
oris(:,2)=rt(oris(:,1));
ori_seg=nan(size(meta,1),1);
for i=1:size(oris,1)
    if i==1
        sell=1:oris(i,1);
    else
        [val,pos]=max(rt(oris(i-1,1):oris(i,1)));
        pos=oris(i-1,1)+pos-1;
        if val>(rt(oris(i,1))+1)
            sell=max(find(meta(:,4)==i,1),pos):oris(i,1);
        else
            sell=[];
        end
    end
    if i==size(oris,1)
        selr=oris(i,1):size(meta,1);
    else
        [val,pos]=max(rt(oris(i,1):oris(i+1,1)));
        pos=oris(i,1)+pos-1;
        if val>(rt(oris(i,1))+1)
            selr=oris(i,1):min(find(meta(:,4)==i,1,'last'),pos);
        else
            selr=[];
        end
    end
    sels=unique([sell,selr]);
    ori_seg(sels,1)=i;
end
dis=nan(size(meta,1),2);
for i=1:size(oris,1)    
    sels=ori_seg==i & ~isnan(rt);
    dis(sels,:)=[meta(sels,3) rt(sels)-oris(i,2)];
    sels=sels& abs(dis(:,1))<=max(ulim) & abs(dis(:,1))>=llim;
    if sum(sels)>4
        oris(i,3)=corr(abs(dis(sels,1)),dis(sels,2),'rows','complete');
        p=robustfit(abs(dis(sels,1)),rt(sels,1));
        oris(i,6)=range(dis(sels,2));
        oris(i,7)=range(dis(sels,1));
    else
        oris(i,[3,6,7])=NaN;
        p=[NaN NaN];
    end
    oris(i,4)=p(2);
end
good_oris=find(meta(oris(:,1),2)<=ip.Results.rtth&oris(:,3)>=0.85);
sels=all(~isnan(dis),2) & ismember(ori_seg,good_oris);
binsize=diff(meta(1:2,1));
if numel(ulim)==1
    sel_fit=sels & abs(dis(:,1))<=ulim&abs(dis(:,1))>=llim;
    medy=accumarray(1+floor((abs(dis(sel_fit,1))-llim)/binsize),dis(sel_fit,2),[],@median);
    medx=accumarray(1+floor((abs(dis(sel_fit,1))-llim)/binsize),abs(dis(sel_fit,1)),[],@median);
    rfit=fitlm(abs(dis(sel_fit,1)),dis(sel_fit,2),'RobustOpts','on');
    plim=rfit.Coefficients.Estimate;
    pmed=robustfit(medx,medy);
    if numel(ip.Results.axis)>0
        axis(ip.Results.axis);
        hold off
        out.s=dscatter(abs(dis(sels,1))/1000,(dis(sels,2)));
        out.s.Children.DisplayName='Loci around ORI';
        hold on;
        %out.sm=scatter(medx,medy,'ko','filled');
        %out.pm=plot([llim ulim],[llim ulim]*pmed(2)+pmed(1),'k--','LineWidth',2,'DisplayName',sprintf('%.2f',1/(pmed(2)*1000))) ;
        out.pl=plot([llim ulim]/1000,[llim ulim]*plim(2)+plim(1),'k-','LineWidth',2,'DisplayName',sprintf('fit%.2f',1/(plim(2)*1000)));
        %out.pl=plot(xlim()/1000,xlim()*plim(2)+plim(1),'k--','LineWidth',2,'DisplayName',sprintf('%.2f',1/(plim(2)*1000)));
        xlabel('distance from next ORI(kb)')
        ylabel('Delta RepTime(loc-ori)')
        title(sprintf('%d/%d',numel(good_oris),size(oris,1)))
        legend('show','Location','best')
    end
else
    axis(ip.Results.axis);
    hold off
    out.s=dscatter(abs(dis(sels,1))/1000,(dis(sels,2)));   
    out.s.Children.DisplayName='Loci around ORI';
    hold on;
    for ui=1:numel(ulim)
        sel_fit=sels & abs(dis(:,1))<=ulim(ui)&abs(dis(:,1))>=llim;
        plim=robustfit(abs(dis(sel_fit,1)),dis(sel_fit,2));
        if numel(ip.Results.axis)>0                                            
            out.pl=plot([llim ulim]/1000,[llim ulim]*plim(2)+plim(1),'-','LineWidth',1,'DisplayName',sprintf('fit%.2f',1/(plim(2)*1000)));
            title(sprintf('%d/%d',numel(good_oris),size(oris,1)))
            legend('show','Location','best')
        end
    end
    xlabel('distance from next ORI(kb)')
    ylabel('\Delta RepTime(loc-ori)')
    xlim([0 20])
    pmed=0;
end
oris(good_oris,end+1)=1;
out.ori_seg=ori_seg;
out.oris=oris;
out.plim=plim;
out.pmed=pmed;
out.rfit=rfit;