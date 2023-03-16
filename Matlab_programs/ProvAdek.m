function mk = ProvAdek(m)
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
if (abs(m)>1.1)  
    m=0;  
end; 
if (abs(m)<0.23)
    m=0;  
end;
mk=real(abs(m));
end