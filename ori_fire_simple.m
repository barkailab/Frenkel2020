function out=ori_fire_simple(rt,meta,varargin)
ip=inputParser;
ip.addParameter('axis',[]);
ip.parse(varargin{:});
%% find oris
sels=~isnan(meta(:,4));
min_dis=accumarray(meta(sels,4),abs(meta(sels,3)),[],@min);
oris=find(abs(meta(sels,3))==min_dis(meta(sels,4)));
sels=~isnan(rt);
[~,x]=findpeaks(rt(sels),find(sels),'MinPeakProminence',1);
replicon_end=false(size(rt));
replicon_end(x)=1;
replicon_end(1)=1;
rep_id=cumsum(replicon_end);
clear replicon_end x
c=0
oris(:,2)=NaN;
for i=unique(rep_id(oris(:,1)))'
    sel_oris=find(rep_id(oris(:,1))==i & rt(oris(:,1))>0);
    if numel(sel_oris)>0
        [~,idx]=min(rt(oris(sel_oris,1)));
        sel_oris=sel_oris(idx);
        sel_i=oris(sel_oris,1)+[-20:20];
        sel_i=sel_i(sel_i>0 & sel_i<size(rt,1));
        oris(sel_oris,2)=corr(rt(sel_i),abs(sel_i-oris(sel_oris,1))','rows','pairwise');
    end
end
selo=oris(:,2)>=.75;
[~,stOri]=mink(meta(oris(:,1),2),5);

rfit=fitlm(meta(oris(selo,1),2),rt(oris(selo,1),1)-median(rt(oris(stOri,1),1)),'RobustOpts','on');
pfi=rfit.Coefficients.Estimate;
out.oris=oris;
out.pfi=pfi;
out.rfit=rfit;
if numel(ip.Results.axis)>0
    axes(ip.Results.axis)
    out.s=scatter(meta(oris(selo,1),2),rt(oris(selo,1),1)-median(rt(oris(stOri,1),1)),'.','displayname',sprintf('ORIs %d',sum(selo)))
    xlabel('RepScore (Yabuki et al. 2002)')
    ylabel('Delta RepTime (ori_{i}-ori_{1st})')
    hold on
    out.l=plot([17 37],[17 37]*pfi(2)+pfi(1),'-','Linewidth',2,'DisplayName',sprintf('fire^{-1}:%.2f',pfi(2)),'Color',out.s.CData)
    legend('show','Location','best')
end
end