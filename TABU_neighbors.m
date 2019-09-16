function [candidate] = neighbors(initial_sol,nt,M,d)
real_value =[];
for integer = -sqrt(M)+1 : 2 : sqrt(M)-1;
    real_value =  [real_value integer];
end
real_value = real_value * d;

for integer =1 : 2*nt
    
    [i,j] = find(real_value==initial_sol(integer,1));
    
    
    if j==1
        temp1 = initial_sol;
        temp1(integer,1) = real_value(2);
        temp2 = initial_sol;
        temp2(integer,1) = real_value(3);
        
        candidate(:,2*integer-1) = temp1;
        candidate(:,2*integer)   = temp2;
        
    elseif j==sqrt(M)
        
        temp1 = initial_sol;
        temp1(integer,1) = real_value(sqrt(M)-1);
        temp2 = initial_sol;
        temp2(integer,1) = real_value(sqrt(M)-2);
        
        candidate(:,2*integer-1) = temp1;
        candidate(:,2*integer)   = temp2;
        
    else
        
        temp1 = initial_sol;
        temp1(integer,1) = real_value(j-1);
        temp2 = initial_sol;
        temp2(integer,1) = real_value(j+1);
        
        candidate(:,2*integer-1) = temp1;
        candidate(:,2*integer)   = temp2;
    end
end
end


