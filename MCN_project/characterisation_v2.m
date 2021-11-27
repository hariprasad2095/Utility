%how many points will the UE cover during characterisation    
clear all;
close all;
x_ref = [1:0.5:10];


%%fileID = fopen('log.txt','w');

x = 37*x_ref;
y = 21*x_ref;
z = ones(1,length(x));
z = 12*x_ref;

RSRP=[];
Throughput=[];

    for i=1:1:length(x_ref)
        %%fprintf(fileID,'%f\n',x_ref);
        [RSRQ_UE(i), RSRP_UE(i)]=BeamManagement(x(i),y(i),z(i));
        Throughput_UE(i)=TPCalculator(RSRQ_UE(i));
        
        RSRP=[RSRP RSRP_UE(i)];
        Throughput=[Throughput Throughput_UE(i)];
        %%fprintf(fileId,'throughput at %d iteration is %f\n',i,Throughput_UE(i));
        %%fprintf(fileId,'RSRP at %d iteration is %f\n',i,RSRP_UE(i));
    end 
    
    tolerable_Throughput=0.9*Throughput(1);
    
    plot(Throughput,'r'); hold on; plot(RSRP,'b');
    
    %finding the interval of tolerable_Throughput 
    for i=2:1:length(x_ref)
        if(tolerable_Throughput>Throughput(i))
            max_int=i;
            min_int=i-1;
            %%fprintf(fileId,'Max_int at %d iteration is %d\n',i,max_int);
            %%fprintf(fileId,'min_int at %d iteration is %d\n',i,min_int);
            break;
        end 
    end 
    
    %fitting a line between these 2 values 
    slope=Throughput(max_int)-Throughput(min_int);
    %%fprintf(fileId,'slope %f\n',slope);

    dist_tol_throughput =(tolerable_Throughput-Throughput(max_int))/slope + max_int;
    %%fprintf(fileId,'dist_tol_throughput  %f\n',dist_tol_throughput);
    %%fclose(fileId);
    %finding RSRP at this dist point 
    y=interp1(RSRP,dist_tol_throughput);


