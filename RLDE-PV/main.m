clear all
close all

rehash

format long g

warning ('off');

run_times = 30;

tic()
for type=1:5
    if type==1
        Max_NFE = 30000;
    elseif type==2
        Max_NFE = 30000;
    elseif type==3
        Max_NFE = 30000;
    elseif type==4
        Max_NFE = 30000;
    elseif type==5
        Max_NFE = 30000;
    end
% for type=[29]
%     if type==20
%         Max_NFE = 30000;
%     elseif type==21
%         Max_NFE = 30000;
%     elseif type==22
%         Max_NFE = 30000;
%     elseif type==27
%         Max_NFE = 30000;
%     elseif type==28
%         Max_NFE = 30000;
%     elseif type==29
%         Max_NFE = 30000;
%     end
    SS=[];
    for i=1:run_times
        [solution ,value ,nfes] = RLDE(@ModelFunction,type,Max_NFE);
        %  [solution ,value ,nfes] = RLDE(@PV_TestFunction_29_s,type,Max_NFE);
        
        SS(i,:) = solution;
        meanss(type,i) = value;
    end
    
    [~,index] = min(meanss(type,:));
    
    [bestvalue,bestindex]=min(meanss(type,:));
    bestsolution = SS(bestindex,:);
    pp(type,:) = [min(meanss(type,:)),max(meanss(type,:)),mean(meanss(type,:)),std(meanss(type,:))];
    
    
    if type==1
        AA1= [SS(index,:),min(meanss(type,:))];
    elseif type==2
        AA2= [SS(index,:),min(meanss(type,:))];
    elseif type==3
        AA3 = [SS(index,:),min(meanss(type,:))];
    elseif type==4
        AA4 = [SS(index,:),min(meanss(type,:))];
    elseif type==5
        AA5 = [SS(index,:),min(meanss(type,:))];
    elseif type==6
        AA6 = [SS(index,:),min(meanss(type,:))];
    end
    
    fprintf('the number of %d function has been finished\n',type);
end
% save AA1
% save AA2
% save AA3
% save AA4
% save AA5
% save AA6
toc()
