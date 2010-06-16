aminoAcidMasses = {}
params = {}
useMonoisotopicMass = False
if useMonoisotopicMass:
   H = aminoAcidMasses['h'] =  1.007825035
   O = aminoAcidMasses['o'] = 15.99491463
   C = aminoAcidMasses['c'] = 12.0000000  
   N = aminoAcidMasses['n'] = 14.0030740  
   P = aminoAcidMasses['p'] = 30.973762   
   S = aminoAcidMasses['s'] = 31.9720707  
else:
   H = aminoAcidMasses['h'] =  1.00794    
   O = aminoAcidMasses['o'] = 15.9994     
   C = aminoAcidMasses['c'] = 12.0107     
   N = aminoAcidMasses['n'] = 14.0067     
   P = aminoAcidMasses['p'] = 30.973761   
   S = aminoAcidMasses['s'] = 32.065     

aminoAcidMasses['G'] = C*2  + H*3  + N   + O 
aminoAcidMasses['A'] = C*3  + H*5  + N   + O 
aminoAcidMasses['S'] = C*3  + H*5  + N   + O*2 
aminoAcidMasses['P'] = C*5  + H*7  + N   + O 
aminoAcidMasses['V'] = C*5  + H*9  + N   + O 
aminoAcidMasses['T'] = C*4  + H*7  + N   + O*2 
aminoAcidMasses['C'] = C*3  + H*5  + N   + O   + S 
aminoAcidMasses['L'] = C*6  + H*11 + N   + O 
aminoAcidMasses['I'] = C*6  + H*11 + N   + O 
aminoAcidMasses['N'] = C*4  + H*6  + N*2 + O*2 
aminoAcidMasses['D'] = C*4  + H*5  + N   + O*3 
aminoAcidMasses['Q'] = C*5  + H*8  + N*2 + O*2 
aminoAcidMasses['K'] = C*6  + H*12 + N*2 + O 
aminoAcidMasses['E'] = C*5  + H*7  + N   + O*3 
aminoAcidMasses['M'] = C*5  + H*9  + N   + O   + S 
aminoAcidMasses['H'] = C*6  + H*7  + N*3 + O 
aminoAcidMasses['F'] = C*9  + H*9  + N   + O 
aminoAcidMasses['R'] = C*6  + H*12 + N*4 + O 
aminoAcidMasses['Y'] = C*9  + H*9  + N   + O*2 
aminoAcidMasses['W'] = C*11 + H*10 + N*2 + O 

aminoAcidMasses['O'] = C*5  + H*12 + N*2 + O*2 
aminoAcidMasses['X'] = aminoAcidMasses['L']  # treat X as L or I for no good reason 
aminoAcidMasses['B'] = (aminoAcidMasses['N'] + aminoAcidMasses['D']) / 2.0  # treat B as average of N and D 
aminoAcidMasses['Z'] = (aminoAcidMasses['Q'] + aminoAcidMasses['E']) / 2.0  # treat Z as average of Q and E 

params['MIN_PEPTIDE_MASS'] =   500.0
params['MAX_PEPTIDE_MASS'] =   6000
params['ALLOWED_MISSED_CLEAVAGES'] = 1
params['MAX_LEN_PEPTIDE']  =   100
params['MIN_LEN_PEPTIDE']  =   6
params['MAX_LEN_DEFINITION'] = 50
params['INITIAL_SEQ_LEN']    = 1000
params['MAX_LEN_PROTEIN']    = 5000
params['PROTON_MASS']        = 1.0072764668

