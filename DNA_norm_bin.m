function [bin_data, bin_meta,out]=DNA_norm_bin(p_all,varargin)
%%is called by parse_full_profile
ip=inputParser;
ip.addParameter('bin_size',500);
ip.addParameter('n_chr',size(p_all,1));
ip.addParameter('norm_nucl',true)
ip.addParameter('calc_ribo_frac',false);
ip.addParameter('calc_mito_frac',false);
ip.addParameter('normVar',false);
ip.addParameter('clean',true);
ip.addParameter('indChr',false);
ip.addParameter('value','mean');
ip.addParameter('Tn5',false);
ip.addParameter('rmvVar',false);
ip.addParameter('cutEnds',false);

ip.addParameter('n_ori',258);


ip.addParameter('format','log2');
ip.parse(varargin{:});
GP=load('group_imp.mat');
n_chr=ip.Results.n_chr;
chr_idx=[0 cumsum(GP.chr_len)];

if ~iscell(p_all) & isnumeric(p_all)
    p_all=mat2cell(p_all,GP.chr_len(1:ip.Results.n_chr),size(p_all,2));
end
if ip.Results.Tn5
    n_vr=size(GP.varRegions,1);
    ty=find(contains(GP.gff.type,'Ty'));
    GP.varRegions{n_vr+1:n_vr+numel(ty),1:3}=[GP.gff.chr(ty),sort([GP.gff.Start(ty) GP.gff.End(ty)],2)+[-200 200]];
    GP.varRegions.name(n_vr+1:end)=GP.gff.type(ty);
    n_vr=size(GP.varRegions,1);
    telo=find(contains(GP.gff.type,'telomere')&contains(GP.gff.dir,'+'));
    GP.varRegions{n_vr+1:n_vr+numel(telo),1:3}=[GP.gff.chr(telo),sort([GP.gff.Start(telo) GP.gff.End(telo)],2)+[-200 0]];
    GP.varRegions.name(n_vr+1:end)=GP.gff.type(telo);  
    n_vr=size(GP.varRegions,1);
    telo=find(contains(GP.gff.type,'telomere')&contains(GP.gff.dir,'-'));
    GP.varRegions{n_vr+1:n_vr+numel(telo),1:3}=[GP.gff.chr(telo),sort([GP.gff.Start(telo) GP.gff.End(telo)],2)+[0 200]];
    GP.varRegions.name(n_vr+1:end)=GP.gff.type(telo); 
end
if ip.Results.cutEnds
    n_vr=size(GP.varRegions,1);
    cutEnd=100000;
    c=0;
    for i=1:16
        c=c+1;
        GP.varRegions{n_vr+c,1:3}=[i 1 cutEnd];
        c=c+1;
        GP.varRegions{n_vr+c,1:3}=[i GP.chr_len(i)-cutEnd+1 GP.chr_len(i)];
    end    
end

bin_data=cell(n_chr,1);
bin_idx=cell(n_chr,1);
total_reads=cellfun(@sum,p_all,'UniformOutput',0);
use_loc=true(sum(GP.chr_len),1);
p_norm=cell(size(p_all));
if ip.Results.norm_nucl & ~ip.Results.indChr
    use_loc(sum(GP.chr_len(1:16))+1:end)=0;
end
if ~ip.Results.normVar
    for i=1:size(GP.varRegions,1)
        use_loc(chr_idx(GP.varRegions{i,1})+[GP.varRegions{i,2}:GP.varRegions{i,3}],1)=false;
    end
end
if ~ip.Results.indChr
    p_temp=cat(1,p_all{1:17});
    base_used=sum(use_loc);
    reads_used=sum(cat(1,total_reads{1:17}))-sum(p_temp(~use_loc,:));
    norm_fac=base_used./reads_used;
    clear p_temp
    for i=1:n_chr
        p_norm{i,1}=p_all{i}.*norm_fac;
    end
    clear norm_fac base_used
else    
    for i=1:n_chr
        use_loc_i=use_loc(chr_idx(i)+1:chr_idx(i+1));
        reads_used(i,:)=sum(p_all{i}(use_loc_i,:));
        fac_i=sum(use_loc_i)./sum(p_all{i}(use_loc_i,:));
        p_norm{i,1}=p_all{i}.*fac_i;        
    end
    clear use_loc_i fac_i
    reads_used=sum(reads_used,1);
end
if ip.Results.rmvVar
    for i=1:16 % we do not remove mitochondrial reads for now
        use_loc_i=use_loc(chr_idx(i)+1:chr_idx(i+1));
        p_norm{i}(~use_loc_i,:)=NaN;
    end
    clear use_loc_i
end
if ip.Results.calc_mito_frac
    out.mito=total_reads{17}./(reads_used+total_reads{17});
    out.cox1=sum(p_all{17}(14000:20000,:),1)./(reads_used+sum(p_all{17}(14000:20000,:),1));
end
clear p_all reads_used;
   
% check normlaisation

chrmean=cellfun(@(x)mean(x,'omitnan'),p_norm,'UniformOutput',false);
figure
imagesc(cat(1,chrmean{:}),[.95 1.05])
clear chr_mean

bin_size=ip.Results.bin_size;
n_sample=size(p_norm{1},2);
bin_chr=cell(n_chr,1);
bin_var=cell(n_chr,1);
bin_idx=cell(n_chr,1);
bin_data=cell(n_chr,1);
for i=1:n_chr
    n_missing=bin_size-mod(size(p_norm{i},1),bin_size);% round(size(p_norm{i},1)/BIN_SIZE)*BIN_SIZE-size(p_norm{i},1);
    if n_missing<=bin_size*0.5
        temp_i=[p_norm{i};nan(n_missing,size(p_norm{i},2))];
        bin_var{i}=mean(reshape([~use_loc(chr_idx(i)+1:chr_idx(i+1));false(n_missing,1)],bin_size,[]))'>0.1;
    else
        temp_i=p_norm{i}(1:end+n_missing-bin_size,:);
        bin_var{i}=mean(reshape(~use_loc(chr_idx(i)+1:chr_idx(i+1)+n_missing-bin_size),bin_size,[]))'>0.1;
    end
    temp_i=reshape(temp_i,bin_size,[],n_sample);
    if ip.Results.clean
        %bad=squeeze(sum(temp_i==0)>0.25*bin_size);
        bad=shiftdim(mean(temp_i>0,1)<=0.75,1);
    end
    if strcmp(ip.Results.value,'mean')
        temp_i=shiftdim(mean(temp_i,1,'omitnan'),1);
    elseif strcmp(ip.Results.value,'sum')
        temp_i=shiftdim(sum(temp_i,1,'omitnan'),1);
    else        
        temp_i=shiftdim(median(temp_i,1,'omitnan'),1);
    end
    if ip.Results.clean
        temp_i(bad)=nan;
    end
    if strcmp(ip.Results.format,'log2')
        bin_data{i}=log2(temp_i);
    elseif strcmp(ip.Results.format,'abs')
        bin_data{i}=temp_i;
    end
    bin_idx{i}=[chr_idx(i)+bin_size/2:bin_size:chr_idx(i+1)]';
    bin_chr{i}=ones(size(bin_idx{i}))*i;    
end
bin_data=cat(1,bin_data{:});
bin_idx=cat(1,bin_idx{:});
bin_chr=cat(1,bin_chr{:});
bin_var=cat(1,bin_var{:});
%ToR = readtable('/home/labs/barkailab/yuliago/Documents/MATLAB/ToR_raz.csv');
ToR =readtable('torSac3.csv');ToR=ToR(:,1:3);
ToR.idx=chr_idx(ToR.chr)'+ToR.start;
[~,idx]=min(abs(bin_idx-ToR.idx'),[],2);
bin_tor=ToR.score(idx);
if ip.Results.n_ori==258
    load('oris.mat','oridb');
    dists=bin_idx-(chr_idx(oridb(:,1))+mean(oridb(:,2:3),2)');
elseif ip.Results.n_ori==410
    dists=bin_idx-GP.oris410.loc';
else    
    dists=bin_idx-GP.oris.loc';   
end
[~,bin_ori]=min(abs(dists),[],2);
bin_ori_dis=dists(sub2ind(size(dists),[1:numel(bin_ori)]',bin_ori));
bin_ori_dis(bin_ori_dis==bin_size/2)=bin_size/2-1;
bin_meta=[bin_idx,bin_tor,bin_ori_dis,bin_ori,bin_chr,bin_var];
if ip.Results.calc_ribo_frac
    out.ribo=sum(p_norm{12}(451000:470000,:));
end
clear bin_idx bin_tor bin_ori_dis bin_ori dists idx p_norm n_tp temp_i BIN_SIZE ToR i bad n_bin
mito=bin_meta(:,1)>chr_idx(17);
bin_meta(mito,2:4)=NaN;
if ~exist('out','var')
    out=[];
end
end
function ori_seg=define_oriseg(meta,bin_size)
ori=abs(meta(:,3))<bin_size/2;
ori_seg=cumsum(ori);
ori_border=false(size(ori_seg));
for i=1:max(ori_seg)-1
    sel_i=find(ori_seg==i);
    [~,idx]=max(meta(sel_i,2));
    ori_border(sel_i(idx))=true;
end
ori_seg=cumsum(ori_border)+1;
end
    

