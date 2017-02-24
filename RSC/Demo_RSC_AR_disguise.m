
clear all;
clc;

%loading path
tem_fd = cd;
par.d_fd          =   [cd '\database\'];
addpath([cd '\utilities\']);

%seting parameter
par.nClass        =   100;
par.nSample       =   7;
par.ID            =   [];
par.nameDatabase  =   'AR_disguise';

%loading data: This is the first experiment of AR with disguise. For the
%second experiment of AR with disguise, please modify the coding loading
%data.
load([par.d_fd 'AR_database']);
Tr_DAT = []; trls = [];
for ci = 1:100
    Tr_DAT = [Tr_DAT Tr_dataMatrix(:,1+7*(ci-1):4+7*(ci-1))];
    Tr_DAT = [Tr_DAT Tt_dataMatrix(:,1+7*(ci-1):4+7*(ci-1))];
    trls   = [trls repmat(ci,[1 8])];
end
clear Tr_dataMatrix Tr_sampleLabels Tt_dataMatrix Tt_sampleLabels;

load([par.d_fd 'AR_database_Occlusion.mat']);
Tt_DAT_sunglass = []; ttls_sunglass = [];
Tt_DAT_scarf = []; ttls_scarf = [];
for ci = 1:100
    Tt_DAT_sunglass = [Tt_DAT_sunglass Tr_dataMatrix(:,1+6*(ci-1):1+6*(ci-1))];
    Tt_DAT_sunglass = [Tt_DAT_sunglass Tt_dataMatrix(:,1+6*(ci-1):1+6*(ci-1))];
    ttls_sunglass   = [ttls_sunglass repmat(ci,[1 2])];
    
    Tt_DAT_scarf    = [Tt_DAT_scarf Tr_dataMatrix(:,4+6*(ci-1):4+6*(ci-1))];
    Tt_DAT_scarf    = [Tt_DAT_scarf Tt_dataMatrix(:,4+6*(ci-1):4+6*(ci-1))];
    ttls_scarf      = [ttls_scarf repmat(ci,[1 2])];
end
clear Tr_dataMatrix Tr_sampleLabels Tt_dataMatrix Tt_sampleLabels;

% Tt_DAT            =  Tt_DAT_sunglass;
% ttls              =  ttls_sunglass;
Tt_DAT              =  Tt_DAT_scarf;
ttls                =  ttls_scarf;

for i = 1:size(Tr_DAT,2)
    tem = reshape(Tr_DAT(:,i),[165 120]);
    tem1 = uint8(imresize(tem,[42 30]));
    O_Tr_DAT(:,i) = tem1(:);
end

Tem_DAT = [];
for i = 1:size(Tt_DAT,2)
    tem = reshape(Tt_DAT(:,i),[165 120]);
    tem1 = uint8(imresize(tem,[42 30]));
    O_Tt_DAT(:,i) = tem1(:);
end

O_Tr_DAT = double(O_Tr_DAT);
O_Tt_DAT = double(O_Tt_DAT);

ID  =  [];
for indTest = 1:size(O_Tt_DAT,2)
    fprintf(['Totalnum:' num2str(size(O_Tt_DAT,2)) 'Nowprocess:' num2str(indTest) '\n']);
    [id] = RSC (O_Tr_DAT,trls,O_Tt_DAT(:,indTest));
    ID   = [ID id];
end

reco = sum(ID==ttls)/length(ttls);
