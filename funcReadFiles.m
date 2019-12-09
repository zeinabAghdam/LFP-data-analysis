function [ str ] = funcReadFiles(  number , dataList  )

% reading the rat subjcet number
if dataList ==1 
 patients = {'contra_06_06_2013_06' ,'contra_06_06_2013_07', 'contra_06_06_2013_08',...
    'contra_06_06_2013_09','contra_25_04_2013_06','contra_25_04_2013_07',...
    'contra_25_04_2013_08','contra_25_04_2013_09','contra_25_04_2013_10'};
else
    patients = {'ipsi_06_06_2013_02','ipsi_06_06_2013_03','ipsi_06_06_2013_04',....
    'ipsi_06_06_2013_05','ipsi_25_04_2013_02','ipsi_25_04_2013_03',...
    'ipsi_25_04_2013_04','ipsi_25_04_2013_05'};
end

if (number>length(patients))||(number<0)
   str = '';
elseif (number==0)
    str = patients;
else
    str = patients(number);
end


end

