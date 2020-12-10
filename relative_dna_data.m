function out=relative_dna_data(data,meta,tp,varargin)
ip=inputParser;
ip.addParameter('ref_data',nan);
ip.addParameter('n_clu',10);
ip.addParameter('cluster',false)
ip.addParameter('adj_rdna',false)
ip.addParameter('sepRef',false)
ip.addParameter('keep',{})
ip.addParameter('strain_id',[])
ip.addParameter('gf2',false)
ip.addParameter('polyfit',true)


ip.parse(varargin{:});
GP=load('group_imp.mat');
chr_idx=[0 cumsum(GP.chr_len)];
% REMOVE VARIABLE REGIONS (E.G. RDNA)
if size(meta,2)>5
    rDNA=meta(:,6)>0;
else
    ty=find(contains(GP.gff.type,'Ty'));
    GP.varRegions{end+1:end+numel(ty),1:3}=[GP.gff.chr(ty),sort([GP.gff.Start(ty) GP.gff.End(ty)],2)+[-500 500]];
    GP.varRegions.name(end-numel(ty)+1:end)=GP.gff.type(ty);
    telo=find(contains(GP.gff.type,'telomere'));
    GP.varRegions{end+1:end+numel(telo),1:3}=[GP.gff.chr(telo),sort([GP.gff.Start(telo) GP.gff.End(telo)],2)];
    GP.varRegions.name(end-numel(telo)+1:end)=GP.gff.type(telo);   
    GP.varRegions{:,5}=mean(GP.varRegions{:,2:3},2)+chr_idx(GP.varRegions{:,1})';
    rDNA=any(abs(meta(:,1)-GP.varRegions{:,5}')<(range(GP.varRegions{:,2:3},2)'/2+1000),2);
end
mito=meta(:,5)==17;
good_sample=mean(isnan(data))<0.1;
gb=all(~isnan(data(:,good_sample)),2);

if numel(ip.Results)>0
    strain_id=ip.Results.strain_id;
else
    strain_id=ones(size(tp));
end
rep_cr=corr(data(gb&~rDNA,:),meta(gb&~rDNA,2),'rows','pairwise')';

if all(isnan(ip.Results.ref_data))
    ref_data=nan(size(data,1),max(strain_id));
    for i=unique(strain_id)
        ssg_i=strain_id==i & ismember(tp,[-5 -1 0]) &good_sample;
        sum(ssg_i)
        if sum(ssg_i)==0
            ssg_i=strain_id==i & (rep_cr>quantile(rep_cr(strain_id==i),1-2/sum(strain_id==i)));
        end
        ref_data(:,i)=mean(data(:,ssg_i),2,'omitnan');
    end   
else
    ref_data=ip.Results.ref_data;
end

if ip.Results.adj_rdna
    good_ref=corr(ref_data,meta(:,2),'rows','pairwise')>-0.1;
    dna_adj=median(ref_data(gb & ~mito &~rDNA,good_ref),'omitnan'); %adj. for ribosomal reads
    ref_data=ref_data-dna_adj;
end
dpl_bin=(ref_data-median(ref_data,2,'omitnan'))>0.8;
if ~ip.Results.sepRef
    ref_data=repmat(median(ref_data,2,'omitnan'),1,max(strain_id));
end    
%%check chromsomes if different strains are available
p=nan(numel(tp),2);
data_rel=nan(size(data));
for i=1:numel(tp)
    ref_data_i=ref_data(:,strain_id(i));
    fit_bin=~isnan(ref_data_i)&~isnan(data(:,i))&~mito&~rDNA & ~dpl_bin(:,strain_id(i)) & meta(:,6)<1;
    if ip.Results.polyfit
        p(i,:)=polyfit(ref_data_i(fit_bin),data(fit_bin,i),1);
    else
        p(i,:)=[1,median(data(fit_bin,i),'omitnan')-median(ref_data_i(fit_bin),'omitnan')];
    end
   	use_bin=~isnan(ref_data_i)&~isnan(data(:,i))&~mito&~rDNA;
    data_rel(use_bin,i)=data(use_bin,i)-ref_data_i(use_bin)*p(i,1)-p(i,2);
    data_rel(mito,i)=data(mito,i)-median(data(use_bin,i)-data_rel(use_bin,i));
    %data_rel(rDNA,i)=data(rDNA,i)-median(data(use_bin,i)-data_rel(use_bin,i));
end
clear fit_bin use_bin ssg_i ref_tp
%rel2=data-ref_data;
data_gs=nan(size(data));
for i=1:max(meta(:,5))
    sel_i=meta(:,5)==i;
    data_gs(sel_i,:)=smoothdata(data_rel(sel_i,:),1,'movmean',11,'omitnan');
end
%% smoothing using meta data
if ip.Results.gf2
    gs2=cell(1,numel(tp));
    gf2=cell(1,numel(tp));
    parfor i=1:size(data_rel,2)
        temp_i=data_rel(:,i);
        temp_i(mito|rDNA)=NaN;
        gs2{i}=nan(size(temp_i));
        gf2{i}=false(size(temp_i));
        for j=find(conv2([1;1;1],1,~isnan(temp_i),'same')>0)'
            sels=abs((1:size(data,1))'-j)<=5 & meta(:,5)==meta(j,5) & ~isnan(meta(:,2)) & ~isnan(temp_i);
            if sum(sels)>7 & range(meta(sels,2))>0.1 & ~isnan(meta(j,2))
                p=robustfit(meta(sels,2),temp_i(sels));
                if numel(lastwarn)==0
                    gs2{i}(j)=meta(j,2)*p(2)+p(1);
                else
                    gs2{i}(j)=mean(temp_i(sels),'omitnan');
                    lastwarn('')
                    gf2{i}(j)=true;
                end
            elseif sum(sels)>3 & (range(meta(sels,2))<0.1 | isnan(meta(j,2)))
                gs2{i}(j)=mean(temp_i(sels),'omitnan');
                gf2{i}(j)=true;
            end
        end
    end
    gs2=cat(2,gs2{:});
    gf2=cat(2,gf2{:});
else
    gs2=[];
end

if numel(ip.Results.strain_id)>0
    long_tc=find(accumarray(strain_id',1)>10);
    % time averageing only for tight time courses use time averages for
    % replication clusters
    data_ts=nan(size(data_gs));%cell(1,numel(long_tc)) %temporal smoothing only for timecourses
    for i=1:numel(long_tc)
        sel_i=strain_id==long_tc(i) & tp>=-5;
        temp=data_gs(:,sel_i);
        temp=temp(:,[1 1:end end]);
        data_ts(:,sel_i)=filter2([.25 .5 .25],temp,'valid');
    end
else
    data_ts=[];
end
clear sel_i;
use_bin=~rDNA & ~mito;
rep_clu=zeros(size(data_gs,1),1);
if ip.Results.cluster
    n_clu=ip.Results.n_clu;
    if numel(long_tc)>0 & ip.Results.cluster>1        
        sel_samples=mean(isnan(data_ts))<0.05;
        rep_clu(use_bin)=kmeans(data_ts(use_bin,sel_samples),n_clu,'MaxIter',500,'Replicates',100,'Options',statset('UseParallel',true));
    else
        sel_samples=mean(~isnan(data_gs(use_bin,:)))>.99;
        rep_clu(use_bin)=kmeans(data_gs(use_bin,sel_samples),n_clu,'MaxIter',500,'Replicates',100,'Options',statset('UseParallel',true));
    end
    [~,rep_order]=sort(accumarray(rep_clu(rep_clu>0),meta(rep_clu>0,2),[],@(x) median(x,'omitnan')));
    rep_clu=changem(rep_clu,[1:numel(rep_order)],rep_order);
    rep_clu(mito)=n_clu+1;
    rep_clu(rDNA)=n_clu+2;
end
out.rel=data_rel;
out.gs=data_gs;
out.ts=data_ts;
out.rep_clu=rep_clu;
out.gs2=gs2;
clear data_ts opts
end