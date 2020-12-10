function [p_af,names]=import_DNA_by_Yoav(path_dir,ending,varargin)
ip=inputParser;
ip.addParameter('outfile','temp.mat');
ip.addParameter('multi',false);
ip.addParameter('subfolder',false)
ip.parse(varargin{:});
GP = load('group_imp.mat');
% LETS find the files directly
% if ~exist('dir_names','var')
%     dir_names=struct2table(dir(path_dir));
%     path_dir=dir_names.folder{1};
%     dir_names=dir_names.name(dir_names.isdir & ~strncmp(dir_names.name,'.',1));
% end
samples=dir([path_dir '/' ending]);

if ip.Results.multi
    [names,~,sample_id]=unique(extractBefore(extractfield(samples,'name')','_R1_'));
    %names=unique(sample_id);
    %[~,sample_id]=ismember(sample_id,names);    
else
    sample_id=1:numel(samples);
    if ~ip.Results.subfolder
        names=extractfield(samples,'name');
    else
        names=extractAfter(extractfield(samples,'folder'),'Sample_');
    end
end
sample_order=table(extractfield(samples,'name')',reshape(sample_id,[],1),reshape(names(sample_id),[],1))
p_all=cell(1,max(sample_id));
tot_reads=cell(1,max(sample_id));
n_chr=17;
GP.chr_len(17)=85779;
sp_reads=cell(1,max(sample_id));
parfor i=1:max(sample_id)
    p_sum=cell(n_chr,1);
    ac=nan(n_chr,601,21);
    files=find(sample_id==i);
    all_dna=zeros(24314210,numel(files));
    sp_reads{i}=nan(2,numel(files))
    for r=1:numel(files)
        fid_i=fopen([samples(files(r)).folder '/' samples(files(r)).name]);
        all_dna(:,r)=fscanf(fid_i,'%g\t');
        fclose(fid_i);
        if exist([samples(files(r)).folder '/' samples(files(r)).name(1:end-4),'_sp.out'],'file')
            sp_reads{i}(:,r)=getpar([samples(files(r)).folder '/' samples(files(r)).name(1:end-4),'_sp.out'])
        end
    end
    sp_reads{i}=sum(sp_reads{i},2);
    tot_reads{i}=sum(all_dna,1);
    all_dna=sum(all_dna,2);
    c=0;
    p=cell(n_chr,2);
    for chr=1:n_chr
        p{chr,1}(:,1)=all_dna(c+1:c+GP.chr_len(chr));
        c=GP.chr_len(chr)+c;
        p{chr,2}(:,1)=all_dna(c+1:c+GP.chr_len(chr));
        c=GP.chr_len(chr)+c;
    end
    for chr=[4,7,13,15,16]
        ac(chr,:,i)=xcorr(movmean(p{chr,1},50),movmean(p{chr,2},50),300);
        ac(chr,:,i)=ac(chr,:,i)/max(ac(chr,:,i));
        %[~,lengths(i,chr)]=findpeaks(ac{i}(chr,301:-1:1),'SortStr','descend','NPeaks',1);
    end
    [~,temp_ext]=findpeaks(mean(ac(:,301:-1:1,i),'omitnan'),'SortStr','descend','Npeaks',1);
    ext_len(i)=max(100,min([250,temp_ext]));
    ext_p=lengthen_profile(p,ext_len(i));
    for chr=1:n_chr
        p_sum{chr,1}=ext_p{chr,1}+ext_p{chr,2};
    end
    p_all(i)={p_sum};
end
p_all=cat(2,p_all{:});
p_final=cell(n_chr,1);
for i=1:n_chr
    for j=1:size(p_all,2)
        p_final{i}(:,j)=p_all{i,j};
    end
end
p_all=p_final;
clear p_final
sp_reads=cat(2,sp_reads{:});
if isfield(ip.Results,'outfile')
    save(ip.Results.outfile,'p_all','names','ext_len','tot_reads','n_chr','sample_order','sp_reads','-v7.3');
end
end

    
% parfor di=1:length(dir_names)
%     try
%         cur_dir_name = dir_names{di};
%         fn = [path_dir '/' cur_dir_name '/mapR1.out'];
%         names{di} = cur_dir_name(length('Sample_X'):end);
%         names{di}(names{di} == '_') = ' ';
%         fprintf('Loading: %s\n%s\n',fn,names{di});        
%         p = sum_profiles_from_files(fn, GP);
%         ext_p = lengthen_profile(p, extend_len);
%         for c=1:16
%             p_cell{di}{c} = ...
%                 ext_p{1,c}' + ext_p{2,c}';
%         end
%     catch
%         names{di} = cur_dir_name(length('Sample_X'):end);
%         names{di}(names{di} == '_') = ' ';
%         
%         for c=1:16
%             p_cell{di}{c} = nan(GP.chr_len(c),1);
%         end
%     end
%     
% end
% p_af=cell(16,1);
% p_cell=cat(1,p_cell{:})';
% for i=1:size(p_cell,1)
%     p_af{i}=cat(2,p_cell{i,:});
% end
% 
% 
% s = zeros(1, length(dir_names));
% for c=1:16
%     s = s+nansum(p_af{c});
% end
% 
% for c=1:16
%     p_af{c} = (p_af{c} ./ s) * sum(GP.chr_len);
% end


function profile =  sum_profiles_from_files(file_name, GP)
profile = cell(2,16);
chr_names=[    
    {'ref|NC_001133|'}
    {'ref|NC_001134|'}
    {'ref|NC_001135|'}
    {'ref|NC_001136|'}
    {'ref|NC_001137|'}
    {'ref|NC_001138|'}
    {'ref|NC_001139|'}
    {'ref|NC_001140|'}
    {'ref|NC_001141|'}
    {'ref|NC_001142|'}
    {'ref|NC_001143|'}
    {'ref|NC_001144|'}
    {'ref|NC_001145|'}
    {'ref|NC_001146|'}
    {'ref|NC_001147|'}
    {'ref|NC_001148|'}
    {'ref|NC_001224|'}];
dirs={'+','-'};
for i=1:16
    for j=1:2
        profile{j,i} = zeros(1,GP.chr_len(i));
    end
end
reads=readtable(file_name,'FileType','text');
for c=1:16
    for j=1:2
        %line = fgetl(fid);
        %C = textscan(line, '%d');
        %C = double(C{1})';
        %fprintf('[Open chip-seq] Chr %d Strand %d, #reads = %d\n',c,(-1)^(j+1),sum(C));
        pos=strcmp(reads.Var1,dirs{j})&strcmp(reads.Var2,chr_names{c});    
        profile{j,c} = accumarray(reads.Var3(pos)+1,1,[GP.chr_len(c) 1])';
    end
end
%fclose(fid);

end

function res_p = lengthen_profile(p, d)
res_p = cell(size(p));
for c=1:size(res_p,1)
    res_p{c,1} = filter2(ones(d,1), p{c,1},'full');
    res_p{c,1} = res_p{c,1}(1:(end-d+1));
    res_p{c,2} = filter2(ones(d,1), p{c,2},'full');
    res_p{c,2} = res_p{c,2}(d:end);
end
end
function rn=getpar(filename)
sp_sam=readtable(filename,'FileType','text','Delimiter','/t','HeaderLines',0,'ReadVariableNames',false);
if ismember('Var1',sp_sam.Properties.VariableNames)
    sp=ismember(sp_sam.Var1,{'NC_001326.1' ,'NC_003421.2','NC_003423.3','NC_003424.3' });
    clear sp_sam
    rn=[sum(sp);sum(~sp)];
else
    rn=[nan;nan]
end
end
