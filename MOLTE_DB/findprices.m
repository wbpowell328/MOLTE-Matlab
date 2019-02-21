function prices = findprices()

    filename = 'PSEGPrices.xlsx'; 
    
%     if (mod(i, 2) == 0) 
%         param = 'B8:B110';
%     elseif (mod(i, 3) == 0) 
%          param = 'B1000:B1100';
%     elseif (mod(i, 5) == 0) 
%          param = 'B10000:B11010';
%     else 
%         param = 'B500:B610';
        
    columnB = xlsread(filename, 'B8:B142400');
    prices = columnB;
    
end 