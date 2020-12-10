xxx = pos*[1000000000;0.5;0.5];
load ../mat' files'/dnaseq_sonication_normalized.mat




n = size(xxx,1);
reads1=reads;
reads1(reads1<0.3)=0.3;
reads1(reads1>2)=2;
data=reads1(:,time_points==-2); %data=reads1(:,time_points==-2);
g1 = mean(reads(:,time_points==-1),2);
data=log2(data./repmat(g1,1,size(data,2)));
async_strains = sample_files(time_points==-2);%can change sample_files with strain_name if you want a shorter name and a mean of strains in the next script.
[dp_num,samples_num] = size(data);
data=(data-repmat(mean(data),dp_num,1))./repmat(std(data),dp_num,1);
xx = x*[1000000000;0.5;0.5];
data2=zeros(n,samples_num);
for i = 1:samples_num
    disp(i);
    %data2(:,i) = smoothData(xx, );
    v=data(:,i);
    v1 = v;
    v1(isnan(v1))=0;
    f1=smooth(v1,10);
    ind1=and(~isnan(v),abs(v1-f1)<1);
    v2 = interp1(xx(ind1),v1(ind1),xxx);
    f2=smooth(xxx,v2,25,'sgolay',2);
    I=f2==0;
    f2(I)=f1(I);
    
    

    f2=csaps(xx,data(:,i),exp(-27.5),xxx);
    data2(:,i)=f2;
end
async_reads = data2;
save ../mat' files'/dnaseq_sonication_async_smooth  pos   async_strains async_reads
