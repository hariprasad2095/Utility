function [UE_X,UE_Y,UE_Z] = Network_UE_Simulator()
    prompt='Base station 1 co-ordinates';
    BS1=[1;2;7];
    bs1_X = BS1(1);
    bs1_Y = BS1(2);
    bs1_Z = BS1(3);
    prompt='Base station 2 co-ordinates';
    BS2=[9;8;9];
    bs2_X = BS2(1);
    bs2_Y = BS2(2);
    bs2_Z = BS2(3);
    prompt='UE_co-ordinates co-ordinates';
    UE=[1,1,0];
    Ue_X = UE(1);
    Ue_Y = UE(2);
    Ue_Z = UE(3);
    
    k=15; % number of steps
    m=1;
    n=k-1;
    UE_X=[];
    UE_Y=[];
    UE_Z=[];
    X_1 = bs1_X;
    Y_1 = bs1_Y;
    X_2 = bs2_X;
    Y_2 = bs2_Y;
    Z_1 = bs1_Z;
    Z_2 = bs2_Z;
    
    while (m <= k-1)&&(n >= 1)
        ue_x = (m*X_2 + n*X_1)/(m+n);
        ue_y = (m*Y_2 + n*Y_1)/(m+n);
        ue_z = (m*Z_2 + n*Z_1)/(m+n);
        UE_X = [UE_X ue_x];
        UE_Y = [UE_Y ue_y];
        UE_Z = [UE_Z ue_z];
        m=m+1;
        n=n-1;
    end
    plot3(UE_X,UE_Y,UE_Z);
end


