%% Load data for both individual iption of ipsi, contra or both%%
function [data, str] = load_data(subjectNumber, opt)
   
   global ldDirec
  
   str = funcReadFiles(  subjectNumber , opt  );
   strFinal = strcat(ldDirec, str{1});
   data = load(strFinal);
end
