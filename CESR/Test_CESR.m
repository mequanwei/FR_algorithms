
function Test_CESR
  fd1=mopen('D:\paper\TSR1\data.bin','rb');
  num=mget(1,'d',fd1);
  ir=mget(1,'d',fd1); 
  ic=mget(1,'d',fd1);
  
  probe60 = mget(ir*ic,'d',fd1); 
  probe60= probe60';

  num = num-1;
  Gallery =[];
  for i = 1:num
    s1 = mget(ir*ic,'d',fd1);
    s1= s1';  
    Gallery=[Gallery s1];
  end  
  
  mclose(fd1);

    
  mprintf('img width: %d; img height: %d \nIteration: ',ic,ir);
  [x,y] = CESR(Gallery,probe60);
  

  disp(find(x>0));

  mprintf('%0.2f ',x(find(x>0)));

  

  [vx,ix]=max(x);
    printf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));
  ax=1:length(x);
  plot2d(ax, x);
endfunction
