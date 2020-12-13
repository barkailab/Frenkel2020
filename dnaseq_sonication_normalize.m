load dnaseq_sonication
ind = sum(reads)>200000;
reads = reads(:,ind);
library_number = library_number(ind);
sample_files=sample_files(ind);


reads=reads*diag(power(mean(reads),-1));%This normalizes all the data to 1 and is essential to not loose everything later where there are filters. 
reads(reads>2)=2;
% reads = normalize(reads,'range',[0,2]);
n=size(reads,2);
strain_name = cell(n,1);
time = cell(n,1);

% To check if a mutant indeed lacks the gene - after running this script,
% look the gene up in sgd and copy its positions  (take a larger range). 
% create a variable with the name of the gene that answers 3 conditions -
% its chromosome is the one you need (x first column), its beginning is
% larger than your lower bound (x second column), and its ending is
% smaller (or equal) to your upper bound (x third column). It will look something like:
% Met7 = and(x(:,1)==15,and(x(:,2)>786000,x(:,3)<800000)); Then,
% imagesc(reads(Met7,:)); At this point you should see a blue spot (x axis
% is samples and y axis is the positions). Now, you should check that the
% samples in the picture are indeed the ones that should lack this gene.
% For this enlarge the picture and use the + tool to find the exact x axis
% positions of your blue bulge. When you have it, take one less and one
% more to test the sample files: sample_files(459:473,:). Check to see that
% you get the same results with strain_name: strain_name(458:474,:). Now,
% all is left is the other side of the same check - make sure that the
% strain you got from your previous search includes only the samples that
% you've found (with the deletions): find(strcmp(strain_name,'met7')) - it
% will print the relevant positions. If you click find(the variable you
% created), you will get the positions of this gene. Then you can test
% whether these positions in x correspond well to the data in Sgd:
% x(19889:19894,:) where the first one is the top result you get, and the
% second one is the bottom. 
% To make a histogram with 100 bins - hist(a1(:),1000) in a hist, if you
% don't want to look at anything higher than x, you write - a1(a1>10)=10;
% hist(a1(:),1000)
% If the reads look all blue, it means that there are outliers sample_filesthat are not
% normalized. So, you have to vizualize: hist(reads(:),1000).

%% fixing names
% The function strrep searches the whole array sample_files for Raz, and in every place 
% changes it to Raz_ - if we didn't call the new variable sample_files also, it would 
% create a new cell array, but here it just updates.
sample_files=strrep(sample_files,'Raz','Raz_');
sample_files=strrep(sample_files,'dp4','dpb4');
sample_files=strrep(sample_files,'Met','Met7');
sample_files=strrep(sample_files,'met','metsic1');
sample_files=strrep(sample_files,'__','_');
sample_files=strrep(sample_files,'with_alfa','-1');
sample_files=strrep(sample_files,'w_alfa','-1');
sample_files=strrep(sample_files,'with_HU','-1');
sample_files=strrep(sample_files,'with_alpha','-1');

sample_files=strrep(sample_files,'before_alfa','-2');
sample_files=strrep(sample_files,'b_alfa','-2');
sample_files=strrep(sample_files,'bef_alpha','-2');
sample_files=strrep(sample_files,'before_alpha','-2');
sample_files=strrep(sample_files,'before_HU','-2');
sample_files=strrep(sample_files,'by4741','wt');
sample_files=strrep(sample_files,'BY4741','wt');
sample_files=strrep(sample_files,'Tail_a','tail-a');
sample_files=strrep(sample_files,'Tail_q','tail-q');
sample_files=strrep(sample_files,'T1_H','Yoav_tos4_HU');
sample_files=strrep(sample_files,'R1_H','Yoav_rad53sml1_HU');
sample_files=strrep(sample_files,'gamma1_H','Yoav_rtt109_HU');
sample_files=strrep(sample_files,'M1_H','Yoav_mec1sml1_HU');
sample_files=strrep(sample_files,'D1_H','Yoav_tos4rtt109_HU');
sample_files=strrep(sample_files,'W1_H','Yoav_sml1_HU');
sample_files=strrep(sample_files,'-INPUT-','_INPUT_');
sample_files=strrep(sample_files,'H3_160_K42Q','h3k42q');
sample_files=strrep(sample_files,'H3_42_K42A','h3k42a');
sample_files=strrep(sample_files,'H3_56_K56A','h3k56a');
sample_files=strrep(sample_files,'Rtt_','Rtt109(nelly)_');
sample_files=strrep(sample_files,'Asf_','Asf1_');
sample_files=strrep(sample_files,'Ada_','Ada1_');
sample_files=strrep(sample_files,'Spt_','Spt7_');
sample_files=strrep(sample_files,'Mec_','Mec1Sml1_');
sample_files=strrep(sample_files,'Sml_','Sml1_');
sample_files=strrep(sample_files,'Rtt109TailA','Rtt109-TailA_-2');
sample_files=strrep(sample_files,'Rtt109TailQ','Rtt109-TailQ_-2');
sample_files=strrep(sample_files,'Sus1tailQ','Sus1-TailQ_-2');
sample_files=strrep(sample_files,'N14-TailR_','N14-TailR');
sample_files=strrep(sample_files,'56QTA','56Q-TailA');
sample_files=strrep(sample_files,'56QTQ','56Q-TailQ');
sample_files=strrep(sample_files,'TQ','TailQ');
sample_files=strrep(sample_files,'Rtt109(nelly)_','Rtt109_');
sample_files=strrep(sample_files,'ada1N11','ada1_N11');
sample_files=strrep(sample_files,'ada2N11','ada2_N11');
sample_files=strrep(sample_files,'ada3N11','ada3_N11');
sample_files=strrep(sample_files,'h3k56a','h3k56a-');




% regexp here returns indexes of letters whithin the string that contain _ here. It can search strings
% for everything.
for i = 1:n
    v=regexp(sample_files{i},'_');
    strain_name{i} = sample_files{i}(1:(v(1))-1);%returns the first word of the string
    if(length(v)>2) %If there are more than 2 _
    time{i} = sample_files{i}((v(end-2)+1):(v(end-1))-1); %This collects the times outomatically, but if they are not written as usual, have to fix. 
    else
        time{i} = '';
    end
end

strain_name(strcmpi(strain_name,'k56a'))={'h3k56a'};
strain_name(strcmpi(strain_name,'k9a'))={'h3k9a'};

%% fixing times
ind = false(n,1);%Here insert libraries that are unsynch data only
ind(library_number==find(strcmp(libraries,'miri_160531_D00257_0229_BC9BLWANXX')))=true;
ind(library_number==find(strcmp(libraries,'nelly_gal_2016')))=true;
ind(library_number==find(strcmp(libraries,'nelly_lowHU_20160721')))=true;
ind(library_number==find(strcmp(libraries,'nelly_Yoav_HU_20160608')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N7_20161009')))=true;
ind(library_number==find(strcmp(libraries,'shai_20170105')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N17_20170208')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N1_20170131')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N19N21N22_20170627')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N24N25N26N27Rnr1_20171008')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N48REDO_20190625/Unsynch')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N36-N45_20180812/Unsynch')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N53-Rtt109TR_20191113')))=true;
ind(library_number==find(strcmp(libraries,'nelly_n53-n54_13112019')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N68_20200803')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N70_20200814')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N71_20200901')))=true;
ind(library_number==find(strcmp(libraries,'nelly_N72_20201020')))=true;


time_points=str2double(time);
%time_points(and(isnan(time_points),ind))=-2;
time_points(ind)=-2;


razyoav_rtt_index = library_number==find(strcmp(libraries,'nelly_Raz_Rtt_20140507'));
yoav_stupid_time_points_system = [39:-3:0,-1];
time_points(razyoav_rtt_index) = yoav_stupid_time_points_system(time_points(razyoav_rtt_index));

razyoav_K9A_K56A_index = library_number==find(strcmp(libraries,'nelly_Raz_ K56AK9A_20150118'));
%yoav_stupid_time_points_system = [39:-3:9,0];
time_points(razyoav_K9A_K56A_index) = [0,33,36,39,9:3:30,0,33,36,39,9:3:30];



%% fixing partial and chromosomal duplications


% The normalizations are done like so: After we run the script, we look for 
% lines (genomic areas) that are yellow. We do this through imagesc(reads); 
% Then, we enlarge the area that interests us and write down the positions in 
% sample_files +-1 to see which library is problematic. Then, we write down X 
% of the patch which intrests us +-1. We then look at these positions
% x(17159:1862) and see which chromosome it is etc. Then, we aply the
% normalization itself. If it is a whole chromosome: 
% ind2=strcmp(strain_name,'ctf18'); ind1=x(:,1)==9;
% factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
% reads(ind1,ind2) = reads(ind1,ind2)*factor;
% If it is a partial duplication:
% ind2=strcmp(strain_name,'h3k56a'); 
% ind1=and(x(:,1)==12,and(x(:,2)>650000 %lower bound%,x(:,3)<817501 %higher bound%));
% factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
% reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'ctf18');
ind1=x(:,1)==9;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n25-sgf29-tailq');%Chromosome 2 duplication
ind1=x(:,1)==2;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n25-tailq-56r');%Chromosome 1 duplication
ind1=x(:,1)==1;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n25-gcn5-56a-taila');%Chromosome 1 duplication
ind1=x(:,1)==1;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n25-gcn5-56a-taila');%Chromosome 16 duplication
ind1=x(:,1)==16;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n27-rtt109-taila');%Chromosome 9 duplication
ind1=x(:,1)==9;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n27-rtt109-taila');%Chromosome 11 duplication
ind1=x(:,1)==11;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n25-56r');% Partial duplication chromosome 4
ind1=and(x(:,1)==4,and(x(:,2)>1102001,x(:,3)<1213000));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n27-tailq-56r');% Partial duplication chromosome 12
ind1=and(x(:,1)==12,and(x(:,2)>657001,x(:,3)<818500));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n26-rtt109-56q');% Partial duplication chromosome 4
ind1=and(x(:,1)==4,and(x(:,2)>884501,x(:,3)<981000));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n25-rtt109-56r-ta');% Partial duplication chromosome 4
ind1=and(x(:,1)==4,and(x(:,2)>519501,x(:,3)<645000));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'n22-rtt109-taila');% Chromosome 1 duplication
ind1=x(:,1)==1;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'rad53-sml1');
ind1=x(:,1)==12;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;


ind2=strcmp(strain_name,'h3k56a'); %partial duplication chr 12
ind1=and(x(:,1)==12,and(x(:,2)>650000,x(:,3)<817501));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'h3k56a'); %partial duplication chr 7
ind1=and(x(:,1)==7,and(x(:,2)>74000,x(:,3)<412001));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'N22-rtt109-56q'); %partial duplication chr 3
ind1=and(x(:,1)==3,and(x(:,2)>138000,x(:,3)<156001));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'N22-rtt109-56q'); %partial deletion chr 3
ind1=and(x(:,1)==3,and(x(:,2)>830000,x(:,3)<840001));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind2=strcmp(strain_name,'N22-tailr'); %Duplication of CIT in chr 4
ind1=and(x(:,1)==14,and(x(:,2)>629000,x(:,3)<630001));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;


selected_list = {'chl4';'h3k42a'};
for i = 1:length(selected_list)
    ind2=strcmp(strain_name,selected_list{i});
    ind1=x(:,1)==1;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;
    
end

selected_list = {'h3k42q';'h3k42a'};
for i = 1:length(selected_list)
    ind2=strcmp(strain_name,selected_list{i});
ind1=x(:,1)==16;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;
    
end

selected_list = {'h3k42q';'h3k42a'};
for i = 1:length(selected_list)
    ind2=strcmp(strain_name,selected_list{i});
ind1=x(:,1)==2;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;
end


selected_list = {'h3k115q';'h3k125q';'h3k18q';'h3k36a';'h3k36q';'h3k37q';'h3k42a';'h3k42q';'h3k64q';'h4k12q';'h4k16a';'h4k16q';'h4k44a';'h3';'h4'};
for i = 1:length(selected_list)
    ind2=strcmpi(strain_name,selected_list{i});
ind1=x(:,1)==14;
factor=(mean(reads(~ind1,ind2)))./(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*diag(factor);
    
end





%ind2=or(strcmp(strain_name,'h3'),strcmp(strain_name,'h4'));
%ind1=x(:,1)==14;
%factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
%reads(ind1,ind2) = reads(ind1,ind2)*factor;

% chr 11 duplication dpb3
ind2=strcmp(strain_name,'dpb3');
ind1=x(:,1)==11;
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

%%
reads(reads>2)=2;
ind=sum(reads<0.7,2)<=(n/2);%instead of 0.5
reads=reads(ind,:);
x=x(ind,:);
ind=sum(reads>1.5,2)<=(n/2);%instead of 0.5
reads=reads(ind,:);
x=x(ind,:);
%%

% chr 9 duplication ctf18
% ind2=strcmp(strain_name,'ctf18');
% ind1=x(:,1)==9;
% factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
% reads(ind1,ind2) = reads(ind1,ind2)*factor;

% partial duplications K56A
ind1 = and(x(:,1)==12,and(x(:,2)>657000,x(:,3)<818001));
ind2=strcmp(strain_name,'K56A');
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

ind1 = and(x(:,1)==7,and(x(:,2)>74000,x(:,3)<412001));
factor=mean(mean(reads(~ind1,ind2)))/mean(mean(reads(ind1,ind2)));
reads(ind1,ind2) = reads(ind1,ind2)*factor;

% remove W303-BY conflicts
%a=(mean(reads(359:364,:)));
%ind1 = a<0.5;
%ind2 = a>=0.5;
%v1=mean(reads(:,ind1),2);
%v2=mean(reads(:,ind2),2);
%ind = abs(v1 - v2)<0.2;
%reads=reads(ind,:);
%x=x(ind,:);

%% In this section we are not normalizing, we are just deleting
% Because there is no signal at all and the area in the genome is very
% small. This (unlike the normalizations) affects all the data
% (we delete the positions) and not just the specific library. These might
% be due to the usage of a different strain from the strain that
% is used in the alignment or duplications. Again, we look at blue lines in
% the matrix reads through imagesc and enlarge the genomic area. We then
% look at X of these positions: x(732000:735002,:) and get the chromosome number etc.
% Then we delete like so: ind=and(x(:,1)==12,and(x(:,2)>732000,x(:,3)< 735002));
% reads=reads(~ind,:);
% x=x(~ind,:); 

% %for the Rif1 library (deletion of Rif1)
% ind=and(x(:,1)==2,and(x(:,2)>751000,x(:,3)< 757501));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% 
% % % for the Mec1 library (deletion of Mec1)
% % ind=and(x(:,1)==2,and(x(:,2)>505500,x(:,3)< 513001));
% % reads=reads(~ind,:);
% % x=x(~ind,:); 
% % 
% % % for the Rad53 library (deletion of Rad53)
% % ind=and(x(:,1)==16,and(x(:,2)>261727,x(:,3)<254192));
% % reads=reads(~ind,:);
% % x=x(~ind,:); 
% 
% %for the Met7 library (deletion of Met7)
% ind=and(x(:,1)==15,and(x(:,2)>787000,x(:,3)< 789001));
% reads=reads(~ind,:);
% x=x(~ind,:); 
% 
% %for the Cac1 library (deletion of RLF2, largest subunit of Cac1)
% ind=and(x(:,1)==16,and(x(:,2)>594001,x(:,3)< 597000));
% reads=reads(~ind,:);
% x=x(~ind,:); 
% 
% 
% %for the ctf18 library (this is ctf18)
% ind=and(x(:,1)==13,and(x(:,2)>422500,x(:,3)< 425001));
% reads=reads(~ind,:);
% x=x(~ind,:); 
% 
% %w3033? met17 or something
% ind=and(x(:,1)==12,and(x(:,2)>732000,x(:,3)< 735002));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %W303 deletion
% ind=and(x(:,1)==7,and(x(:,2)>0,x(:,3)<9001));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% 
% %w303 - deletion
% ind=and(x(:,1)==14,and(x(:,2)>765500,x(:,3)<779001));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %w303 deletion
% ind=and(x(:,1)==8,and(x(:,2)>93000,x(:,3)<95001));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %YCL019W  : Retrotransposon TYA Gag
% ind=and(x(:,1)==3,and(x(:,2)>=84001,x(:,3)<= 93002));
% reads=reads(~ind,:);
% x=x(~ind,:);
%  
% %bar1
% ind=and(x(:,1)==9,and(x(:,2)>= 320501,x(:,3)<= 326001));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% ind=and(x(:,1)==8 ,x(:,2)>= 525501);% 
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% ind=and(x(:,1)==1,and(x(:,2)>=13001,x(:,3)<= 28001));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %ATG1
% ind=and(x(:,1)==15,and(x(:,2)>=159001,x(:,3)<= 162001));
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %Beginning of chromosome 3 Gex1 and Vba3
% ind=and(x(:,1)==3 ,x(:,2)<= 11001); 
% reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %mrc1
% ind=and(x(:,1)==3,and(x(:,2)>=18501,x(:,3)<= 22501));
%  reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %fkh2
% ind=and(x(:,1)==14,and(x(:,2)>= 495501,x(:,3)<= 496501));
%  reads=reads(~ind,:);
% x=x(~ind,:);
% %trp1
% ind=and(x(:,1)==4,and(x(:,2)>= 461501,x(:,3)<= 462501));
%  reads=reads(~ind,:);
% x=x(~ind,:);
% %clb5
% ind=and(x(:,1)==16,and(x(:,2)>= 773501,x(:,3)<= 775501));
%  reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %fkh1
% ind=and(x(:,1)==9,and(x(:,2)>= 100501,x(:,3)<= 102501));
%  reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %ura3
% ind=and(x(:,1)==5,and(x(:,2)>= 116001,x(:,3)<= 117001));
%  reads=reads(~ind,:);
% x=x(~ind,:);
% 
% %N45-56Q-SAGA-TQ1
% %lys2
% ind=and(x(:,1)==2,and(x(:,2)>= 469500,x(:,3)<= 474501));
%  reads=reads(~ind,:);
% x=x(~ind,:);

%%
reads=reads*diag(power(mean(reads),-1));% normalize reads number to mean 1 again after it changed a bit during normalizations
% reads = normalize(reads,'range',[0,2]);
%% Correction of names after libraries already exist
strain_name = lower(strain_name);
strain_name(and(strcmp(libraries(library_number),'data_wt_feb14'),strcmp(strain_name,'wt')))={'wt_yoav'};
strain_name(and(strcmp(libraries(library_number),'miri_141225_D00257_0160_BC5RPDACXX'),strcmp(strain_name,'fkh')))={'fkhsic1'};
strain_name(and(strcmp(libraries(library_number),'miri_141201_SN808_0168_BC5RU4ACXX'),strcmp(strain_name,'by')))={'wt'};
strain_name(and(strcmp(libraries(library_number),'miri_151006_D00257_0196_AC7V1NANXX'),strcmp(strain_name,'fkh')))={'fkh_HU'};
strain_name(and(strcmp(libraries(library_number),'miri_160216_D00257_0214_AC8BK1ANXX'),strcmp(strain_name,'fkh')))={'fkh_240min'};
strain_name(strcmp(strain_name,'tef'))={'tef2sic1'};
strain_name(strcmp(strain_name,'cdc6'))={'metcdc6'};
strain_name(strcmp(strain_name,'rnr'))={'rnr1'};

%[a1,~,strain_id] = unique(strcat(libraries(library_number),'+++',lower(strain_name)));
% a1=regexp(a1,'+++','split');
% n = size(a1,1);
% strain_list = cell(n,1);
% for i = 1:n
%     strain_list{i}=a1{i}{2};
% end
 
%% Re-organization of the data in an ordered way

[strain_list,~,strain_id] = unique(lower(strain_name));

  [~,new_order]=sort(strain_id*1000+time_points);
 strain_id = strain_id(new_order);
 reads = reads(:,new_order);
 sample_files = sample_files(new_order);
 time_points = time_points(new_order);
library_number=library_number(new_order);
strain_name=strain_name(new_order);

save dnaseq_sonication_normalized reads x library_number libraries sample_files time_points strain_list strain_id strain_name

%library_number libraries
