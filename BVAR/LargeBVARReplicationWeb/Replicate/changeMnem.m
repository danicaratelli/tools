function seriesPl = changeMnem(seriesPl);



seriesPl = strrep(seriesPl,'CES002','EMPL');
seriesPl = strrep(seriesPl,'PUNEW','CPI');
seriesPl = strrep(seriesPl,'FYFF','FFR');
seriesPl = strrep(seriesPl,'PSM99Q','COMM PR');
seriesPl = strrep(seriesPl,'FMRRA','TOT RES');
seriesPl = strrep(seriesPl,'FMRNBA','NBORR RES');
seriesPl = strrep(seriesPl,'FM2','M2');
seriesPl = strrep(seriesPl,'A0M051','INCOME');
seriesPl = strrep(seriesPl,'A0M224_R','CONSUM');
seriesPl = strrep(seriesPl,'IPS10','IP');
seriesPl = strrep(seriesPl,'A0m082','CAP UTIL');
seriesPl = strrep(seriesPl,'LHUR','UNEMPL');
seriesPl = strrep(seriesPl,'HSFR','HOUS START');
seriesPl = strrep(seriesPl,'PWFSA','PPI');
seriesPl = strrep(seriesPl,'GMDC','PCE DEFL');
seriesPl = strrep(seriesPl,'CES275','HOUR EARN');
seriesPl = strrep(seriesPl,'FM1','M1');
seriesPl = strrep(seriesPl,'FSPCOM','S&P');
seriesPl = strrep(seriesPl,'FYGT10','TB YIELD');
seriesPl = strrep(seriesPl,'EXRUS','EXR');
