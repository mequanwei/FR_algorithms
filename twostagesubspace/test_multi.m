%test program
tmp=100;
M=tmp;
refine_time=0
for i=1:size(M,2)
  current_cputime=cputime;
  disp(['M=',num2str(M(i))])
  res_correct_rate_3refine(i,:) =TPTSR_multi(A,40,10,M(i),refine_time)  ;
%   save(['M',num2str(M(i)),'.mat'],'res_correct_rate_1refine');
  disp(['CPUTIME=',num2str(cputime-current_cputime)])
end
 save(['ORL correct_rate_3refine','.mat'],'res_correct_rate_3refine');

refine_time=2
for i=1:size(M,2)
  current_cputime=cputime;
  disp(['M=',num2str(M(i))])
  res_correct_rate_2refine(i,:) =TPTSR_multi(A,40,10,M(i),refine_time)  ;
%   save(['M',num2str(M(i)),'.mat'],'res_correct_rate_1refine');
  disp(['CPUTIME=',num2str(cputime-current_cputime)])
end
 save(['ORL correct_rate_2refine','.mat'],'res_correct_rate_2refine');
 
refine_time=1
for i=1:size(M,2)
  current_cputime=cputime;
  disp(['M=',num2str(M(i))])
  res_correct_rate_1refine(i,:) =TPTSR_multi(A,40,10,M(i),refine_time)  ;
%   save(['M',num2str(M(i)),'.mat'],'res_correct_rate_1refine');
  disp(['CPUTIME=',num2str(cputime-current_cputime)])
end
 save(['ORL correct_rate_1refine','.mat'],'res_correct_rate_1refine');
 
refine_time=0
for i=1:size(M,2)
  current_cputime=cputime;
  disp(['M=',num2str(M(i))])
  res_correct_rate_0refine(i,:) =TPTSR_multi(A,40,10,M(i),refine_time)  ;
%   save(['M',num2str(M(i)),'.mat'],'res_correct_rate_1refine');
  disp(['CPUTIME=',num2str(cputime-current_cputime)])
end
 save(['ORL correct_rate_0refine','.mat'],'res_correct_rate_0refine');
