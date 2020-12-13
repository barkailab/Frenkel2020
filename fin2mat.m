function [ ] = fin2mat(  )
%FIN2MAT Summary of this function goes here
%   Detailed explanation goes here

base_path = '/home/labs/barkailab/LAB/data/SEQ';

%DPN - you don't have a reason to touch the next 10 lines
libraries = {'miri_121112_SN808_0102_AD1FKRACXX';...
'miri_121227_SN808_0106_BC0G4HACXX';...
'miri_130311_SN808_0123_BD1FJGACXX';...
'miri_130410_SN808_0126_AD1FN9ACXX';...
'miri_130611_SN808_0139_BC1TYDACXX';...% paired end
'miri_131009_SN808_0151_AD2D13ACXX';...
'miri_131231_SN808_0157_AD2G1AACXX';...
'miri_140205_D00257_0132_BD2H07ACXX'};
saveLib(base_path,libraries,'dnaseq_dpn',false);



%SONICATION - ADD all new data to here!
% Add the new library name to ths cell array
libraries ={'miri_141201_SN808_0168_BC5RU4ACXX';...
'miri_141225_D00257_0160_BC5RPDACXX';...
'miri_150216_SN808_0174_BC5RHMACXX';...
'miri_151006_D00257_0196_AC7V1NANXX';...
'miri_160216_D00257_0214_AC8BK1ANXX';...
'nelly_dpb_20160419';...
'nelly_rnr1_20160322';...
'miri_160531_D00257_0229_BC9BLWANXX';...
'nelly_gal_2016';...
'miri_histone_modifications_4strains_160724';...
'nelly_lowHU_20160721';...
'nelly_Raz_Rtt_20140507';...
'nelly_Yoav_HU_20160608';...
'nelly_Raz_ K56AK9A_20150118';...
'nelly_N5_20161010';...
'nelly_N7_20161009';...
'nelly_N9_20161009';...
'shai_20170105';...
'nelly_N12_20170108';...
'nelly_N17_20170208';...
'nelly_N18_20170208';...
'nelly_N1_20170131';...
'nelly_N14_20170123';...
'nelly_N19N21N22_20170627';...
'nelly_N24N25N26N27Rnr1_20171008';...
'nelly_N39N40_TRTC_20181026';...
'nelly_N48REDO_20190625/Unsynch';
'nelly_N48REDO_20190625/Time courses';...
'nelly_N36-N45_20180812/Time_courses';...
'nelly_N36-N45_20180812/Unsynch';...
'nelly_N53-Rtt109TR_20191113';
'nelly_n53-n54_13112019';
'nelly_N68_20200803'; % paired end
'nelly_N70_20200814'; % paired end
'nelly_N71_20200901'; % paired end
'nelly_N72_20201020'}; % Paired end
% 'nelly_20180307_nextera_calib';...
% 'nelly_20180508_HMTn5_calib';...
% 'nelly_20180531_finalbuffer';...
% 'nelly_eyal_Tn5_SPRI_cleanup_upscaled';...
% 'nelly_Eyal_Tn5';...
% 'nelly_N47test_20190331';
% 'nelly_pairedend_N36WT';};

saveLib(base_path,libraries,'dnaseq_sonication',false);% now its without yoav's data

end

function []  = saveLib(base_path,libraries,matfilename,ADD_YDATA)
reads=[];
sample_files=[];
library_number=[];
for i = 1:length(libraries)
   % disp(libraries{i});
    a=dir(sprintf('%s/%s/Unaligned_fastq/',base_path,libraries{i}));
    if(isempty(a))
         fprintf('could not find files in directory %s\n',libraries{i});
    end
    m=length(a);
    for j = 3:m
        sample = a(j).name;
        b=dir(sprintf('%s/%s/Unaligned_fastq/%s/*.fin',base_path,libraries{i},sample));
        if(isempty(b))
            fprintf('could not find files in sample %s\n',sample);
        end
        sc=struct2cell(b);
        s=sc(1,:)';
        sample_files = [sample_files;s];
        n=length(b);
        for k = 1:n
            fid = fopen(sprintf('%s/%s/Unaligned_fastq/%s/%s',base_path,libraries{i},sample,b(k).name));
            C=textscan(fid,'%d %d %d %d');
            fclose(fid); 
            reads = [reads,C{4}];
            library_number=[library_number;i];
        end
        if(isempty(b))
           fprintf('lib: %s sample: %s\n',libraries{i},sample) ;
        end
    end
end
reads = double(reads);
x=double([C{1},C{2},C{3}]);
ind = x(:,1)<17;
reads = reads(ind,:);
x=x(ind,:);
if(ADD_YDATA)    
    [x,reads,libraries,library_number,sample_files] = addYoavsData(x,reads,libraries,library_number,sample_files);
end

save (matfilename, 'reads','sample_files','x','libraries','library_number');    
end

function[x,reads,libraries,library_number,sample_files] = addYoavsData(x,reads,libraries,library_number,sample_files)

n = size(x,1);
prev_num = size(libraries,1);
yoav_data_path='/home/labs/barkailab/yoavv/Barkai/Projects/HappyCycle2/ArticleFigures/replication/alpha_factor_data_organized';
a=dir(yoav_data_path);
for i = 3:length(a)
    mat_file = a(i).name;
    disp(mat_file);
    libraries = [libraries ;mat_file(1:end-4)];
    A=load(sprintf('%s/%s',yoav_data_path,mat_file));
    B=fieldnames(A);
    if(length(B)>1)
        for j = 1:length(B)
            if(regexp(B{j},'input'))
                C=B{j};
            end
        end
        t=A.all_times;
        A=A.(C);
        m = size(A{1},2);
        Z = zeros(n,m);
        for j = 1:n
            my_chr=x(j,1);
            my_vec=x(j,2):min(x(j,3),length(A{my_chr}));
            
            %fprintf('point %d chr %d pos:%d-%d\n',j,my_chr,my_vec(1),my_vec(end));
            Z(j,:) = nansum(A{my_chr}(my_vec,:));
        end
        samples = cell(m,1);
        for j = 1:m
            samples{j}=sprintf('%s_%d_NNNNNN_L000.fin',mat_file(6:end-10),t(j));
        end
        %fprintf('libnum %d i %d',prev_num,i);
        library_number=[library_number;(prev_num+i-2)*ones(m,1)];
        sample_files = [sample_files;samples];
        reads=[reads,Z];
    end
end


end



