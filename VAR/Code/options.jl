## Author: Daniele Caratelli
#email: danicaratelli@gmail.com (don't hesitate to contact if you are
#confused, at worst I will not answer you)

#description: options for the VAR, this is the only 
#             file you need to modify in order to
#             use the VAR

## ------------------ Options -------------------- ##

   #VAR options
p = 4;                                   #Number of lags
date_format = "mm/d/yyyy";               #Date format
start_date = "01/01/1960";                  #Start date
end_date = "01/01/2009";                    #End date
Vars = ["gdp" "price" "ffr"];            #Variable names (in order)

#Variable names (in order)
Mnems = ["GDPC1","GDPDEF","FEDFUNDS"]; #series to use

Var_names = ["GDP","Inflation","Interest Rate"];

dates_mnems = "date";               #Date mnemonic

Transformations = ["log",
                   "log",
                   "lin];   		#Transformations applied,
                                       #if you want to add 
							#transformations, 
                                       #change the 										#"Data_organizer.jl" accordingly
                   
	#IRF options
ndraws = 1000;              #Number of draws
nperiods = 24;              #Number of shock horizons
conf_int = 66;              #Confidence interval in plot

#Folder names
irf_fold = "IRFs";
fore_fold = "fores";

#Notes: make sure you have three folder at the same level:
#       1) Code (or VarCode), contains all the code
#       2) Data, contains excel sheet with Data
#       3) Outdata, where all is outputted

## ------------------ Options -------------------- ##
