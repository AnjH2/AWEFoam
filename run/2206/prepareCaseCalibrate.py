#function to calculate the watervapour content 
def calVapour(m_KOH,T_KOH,P_KOH):
    import math;
    k1=-0.01508;
    k2=-0.0016788;
    k3=2.25887e-5;
    k4=-0.0012062;
    k5=5.6024e-4;
    k6=-7.8228e-6;
    k7=35.4462;
    k8=-3343.93;
    k9=-10.9;
    k10=0.0041645;
    f=k1*m_KOH+k2*pow(m_KOH,2)+k3*pow(m_KOH,3)+(1+k4*m_KOH+k5*pow(m_KOH,2)+k6*pow(m_KOH,3))*(k7+k8/T_KOH+k9*math.log10(T_KOH)+k10*T_KOH)
    result=(pow(10,f)/P_KOH)
    return result

import numpy as np
import sys
import os
import datetime
import decimal
#########################################################################################################################
#Assing input VARibles
VAR=sys.argv[2] 	#Input when script is called
VAR=VAR.split(',') 	#Splitting input into a list

CASE=sys.argv[1]	#Case location
NPro=int(sys.argv[3])	#Number of cores

Temp=91;		#system temperature
#tauc=1;		#not used anymore	
KPE=763000000.0;	#electrode permeability
I0_Scale=float(VAR[5]); #Value used to scale I0 for both Ne and Pe
VapourCase="true";	#bool to for water vapour inclussion
molality_KOH=7.65;    #value used if Vaporcase true, concentration in molality[mol/kg] and not molarity[mol/m³].
CModel="vogt";		#Bubble coverage model
nuModel="linearMu";	#viscosity model
GG=-9.81#VAR[4]		#Gravity defeault to X direction!


R=8.3145;		#gas constant
## Inputs

#defining domain values
TEL=1.8			#Electrode thickness, currently symmetric[mm]
TDI=0.85		#diaphrame thickness[mm]
WC=30			#Width of cell[mm]
LC=float(VAR[4])	#Length of cell[mm]
#FIELD VALUES
r_p_E	= 0.02e-3	#electrode pore radius [m], maybe used for reynolds number
eps_E	= 0.97#0.95	#electrode porosity
eps_S	= 0.5		#seperator/diaphrame porosity
r_p_S	= 65e-9#130e-9	#seperator/diaphrame pore radius [m], maybe used for reynolds number
tau_s	= float(VAR[3]);#seperator/diaphrame tortuosity, strongly related to ohmic losses


##defining UNe for simulation
UNE_start=-1.4		#Start potential
diff=0.1;		#Potential increase for each step.
nsteps=1+7;		#Number of steps in simulation, more steps longer simulation
stepTime=20; 		#Seconds at constat potential.
ExStartTime=20; 	#Extra startTime.
mstep_ratio=1/4;	#Time ratio for increassing potential
potentialSteps=[];	#List to store potential intervals
timeSteps=[];		#List to store time intervals
for i in range(nsteps): #loop to create potential schedual 
	potentialSteps.append(UNE_start-diff*i);
	timeSteps.append(ExStartTime+(i+1)*stepTime)
	timeSteps.append(timeSteps[i*2]+stepTime*mstep_ratio)
#To be able to compare with experiments some values will be replaced with experimental values.
#replacing with relevant voltages from experiment
base="~/experimental_results/";
expL="exp_L30W30T91Q100.txt"
a,b=np.loadtxt(os.path.expanduser(base+expL), unpack=True)
DU=list(b)
DI=list(a)
for i in DU:
	for j in range(1,len(potentialSteps)):
		diff1=abs((potentialSteps[j-1]+i)/i);
		diff2=abs((potentialSteps[j]+i)/i);
		if (diff1<=diff2):
			potentialSteps[j-1]=-i
			break
#inlet or initial conditions
VolFlow=200 	#inlet flowrate(total) [mL/min]
velocityY=0
velocityX=VolFlow/(1000*1000*60)*1/(30*1.8*2*10e-6)#value from Marcus   #(Vol*Ntime)/(7.2e-5*3600)
velocityZ=0
	
#inlet and initial concentration mol/m³. electrolyte
CH2E=0.001;	#should be low, but not important
CO2E=0.001;	#should be low, but not important
CH2OE=48810;	#can be calculated based on OH and density.
COHE=6800;	#Very important because of conducivity depends on it.
#inlet mol/m³. gas
PNUM=101325; #abselute pressure in cell, model is not ready for changeing value.
if VapourCase:
	p_water=calVapour(molality_KOH,Temp+273.15,PNUM/1e5);#waters partial pressure
else:
	p_water=0;#waters partial pressure
ALPHAIN=(CH2G+CO2G)*R*(Temp+273.15)/(PNUM-p_water*1e5);#inlet and initial volumefraction


#properties for models [CASE/constant/*properties 
	
#H.Vogt
JLim=300e3; #limiting current density for vogt model
	
	
#butlerVolmer
I0NE=7e-1; #exchange current density coefficient not the final value as it will be scaled
I0PE=8e-3; #exchange current density coefficient not the final value as it will be scaled
Ea0_Pe=46156.12; #constant to capture the temperature dependence on I0Pe
Ea0_Ne=40e3;#constant to capture the temperature dependence on I0Pe
aPe=0.61;# transfer coefficient at Pe
aNe=0.715;# transfer coefficient at Ne
	
#Time Settings
#Initial time step
DELTAT        = 0.0005 #time step
ENDTIME       = nsteps*stepTime #end time
WRITEINTERVAL = stepTime*0.2 #writeing interval

MAXCO	= 0.15; 	#max courant number
ADJUSTTIME = 'yes'; 
MAXDT	=	0.05; 	#max timestep
	
#Number of processors
NProcessor=NPro #number processors
print('Number of processors '+str(NProcessor))

#mesh calculations
SS1=2			#Scalling of cells in flow direction, high number less cells

Ny=int(VAR[0])	#Number of cells in the wall normal diretion
SDi=2		#extra cells in diaphrame
t3=TEL*2+TDI	#thickness of cell
#calculating cells in each area.
NyPe=round(TEL/t3*Ny)
NyPe2=round(NyPe/2)
NyDi=round(TDI/t3*Ny*SDi)
NyDi2=round(NyDi/2)

#loop to find growth rate
print("NyDi="+str(NyDi))
print("NyPe="+str(NyPe))
tCW=TEL/(NyPe*3)#smallest wall normal cell
print("smallest wall normal cell tCW="+str(tCW)+"\n")
if tCW>(TDI/NyDi):
	print("\n****\nthe goal of smallest cell is larger than an uniform mesh in the membrane\n changing the goal to uniform membrane mesh\n****\n")
	tCW=TDI/NyDi
error=10
r_guess=1.001
r_Pe=r_guess
i=0
i_max=1e6
while ((error>1e-8)):
	r_guess=(1-(TEL/2)/tCW*(1-r_Pe))**(1./(NyPe2))
	error=abs(r_guess-r_Pe)
	r_Pe=r_guess
	i=i+1
	if(i>i_max):
		print("to many iteration")
		break
error=10
r_Di=r_Pe 
i=0     
while ((error>1e-8)):
	r_guess=(1-(TDI/2)/tCW*(1-r_Di))**(1./(NyDi2))
	error=abs(r_guess-r_Di)
	r_Di=r_guess
	i=i+1
	if(i>i_max):
		print("to many iteration")
		break
r_Pe=r_Pe**(NyPe2-1)
r_Di=r_Di**(NyDi2-1)


#A file to recap the case#
case=CASE

os.system("touch "+case+"/caseInformation.txt")
f=open(case+"/caseInformation.txt","a")	
f.write("caseName:i0Scale:"+str(I0_Scale)+"_Temp[C]:"+str(Temp)+"_Taus:"+str(tau_s)+"_1\n")
f.write("date:"+str(datetime.datetime.now())+"\n")
f.write("caseProperties:")
f.write("cellPotential[V]:"+str(UNE_start)+":"+str(UNE_start-diff*(nsteps-1))+"\n")
f.write("temperature[K]:"+str(Temp+273.15)+"\n")
f.write("KoHConcentration[mol/m³]:"+str(COHE)+"\n")
f.write("inletFlowRate[mL/min]-[m³/s]:"+str((velocityX*(TEL*2+TDI)/1000*WC/1000)*60*1000000)+":"+str((velocityX*(TEL*2+TDI)/1000*WC/1000))+"\n")
f.write("electrode[mm]:"+str(TEL)+":"+str(TEL)+"\n")
f.write("membransize[mm]:"+str(TDI)+"\n")
#f.write("totalNumberOfCells:"+__)
f.write("domainSizeWC[mm]:"+str(LC)+":"+str(TEL+TEL+TDI)+":"+str(WC)+"\n")
#f.write("baseCellSize[mm]:"+__+":"+__+":"+__)
f.write("otherProperties:\n")
f.write("alpha_C_Ne/alpha_A_Pe:"+str(aNe)+"/"+str(aPe)+str("\n"))
f.close()

print ("Current directory is: ")
os.system('pwd')

startFrom="0"
	

#prepareing the case folder	
os.system('rm -rf '+ case+'/0 '+case+'/sets '+case+'/log.*')
os.system('cp -r '+ case+'/0.orig '+ case+'/0')

	
#Control dict#########################################################################################

controlDictOrig=case+"/system/controlDict.orig"
controlDict=case+"/system/controlDict"
os.system('sed -e s/DELTAT/'+str(DELTAT)+'/g  '+controlDictOrig+' > '+controlDict)
os.system('sed -i s/STARTFROM/'+startFrom+'/g  '+controlDict)
os.system('sed -i s/ENDTIME/'+str(ExStartTime+ENDTIME)+'/g  '+controlDict)
os.system('sed -i s/WRITEINTERVAL/'+str(WRITEINTERVAL)+'/g  '+controlDict)

os.system('sed -i s/ADJUSTTIME/'+ADJUSTTIME+'/g  '+controlDict)
os.system('sed -i s/MAXDT/'+str(MAXDT)+'/g  '+controlDict)
os.system('sed -i s/MAXCO/'+str(MAXCO)+'/g  '+controlDict)
#os.system('sed -i s/DEPTH/'+str(ZMAX-ZMIN)+'/g  '+controlDict)
	

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
#insert in blockMesh
	
blockMeshDictOrig=case+"/system/blockMeshDict.orig"
blockMeshDict=case+"/system/blockMeshDict"
	    
os.system('sed -e s/NY/'+str(NZ)+'/g  ' + blockMeshDictOrig+' > '+blockMeshDict)
os.system('sed -i s/SS1/'+str(SS1)+'/g  ' +blockMeshDict)
os.system('sed -i s/SS2/'+str(SS1)+'/g  ' +blockMeshDict)
os.system('sed -i s/TPE/'+str(TEL)+'/g  ' +blockMeshDict)
os.system('sed -i s/TNE/'+str(TEL)+'/g  ' +blockMeshDict)
os.system('sed -i s/TDI/'+str(TDI)+'/g  ' +blockMeshDict)
os.system('sed -i s/LC/'+str(LC)+'/g  ' +blockMeshDict)
os.system('sed -i s/WC/'+str(WC)+'/g  ' +blockMeshDict)
os.system('sed -i s/RDI/'+str(r_Di)+'/g  ' +blockMeshDict)
os.system('sed -i s/RPE/'+str(r_Pe)+'/g  ' +blockMeshDict)
os.system('sed -i s/SDI/'+str(SDi)+'/g  ' +blockMeshDict)
	    
os.system('echo running blockMesh...')
os.system('blockMesh -case '+ case+' > '+case+'/log.blockMesh')
	    
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
#insert into /0	   
 
#velocity initital conditions and BCs
initial=case+"/0/"
os.system('sed -i s/VEX/'      + str(velocityX)     + '/g  ' + initial+'U')
os.system('sed -i s/VEY/'      + str(velocityY)     + '/g  ' + initial+'U')
os.system('sed -i s/VEZ/'      + str(velocityZ)     + '/g  ' + initial+'U')
os.system('sed -i s/TEMP/'	+str(273.5+Temp)+'/g '+initial+'T')

#potential initital conditions and BCs
#nsteps, msteps_ratio, difference and UNE_start is defined above
#potentialSteps	timeSteps
ramptime=ENDTIME/nsteps*mstep_ratio
os.system('sed -i s/UNE0/'+str(potentialSteps[0])+'/g '+initial+'UNe')
os.system('sed -i s/ENDTIME/'+str(timeSteps[0])+'/g '+initial+'UNe')
os.system("sed -i 's/UNE1/"+str(potentialSteps[0])+")(ENDTIME UNE1/g' "+initial+'UNe')
for j in range(1,nsteps-1):
	os.system('sed -i s/ENDTIME/'+str(timeSteps[j*2-1])+'/g '+initial+'UNe')
	os.system("sed -i 's/UNE1/"+str(potentialSteps[j])+")(ENDTIME UNE1/g' "+initial+'UNe')
	os.system('sed -i s/ENDTIME/'+str(timeSteps[j*2])+'/g '+initial+'UNe')
	os.system("sed -i 's/UNE1/"+str(potentialSteps[j])+")(ENDTIME UNE1/g' "+initial+'UNe')
j=j+1;
os.system('sed -i s/ENDTIME/'+str(timeSteps[j*2-1])+'/g '+initial+'UNe')
os.system("sed -i 's/UNE1/"+str(potentialSteps[j])+")(ENDTIME UNE1/g' "+initial+'UNe')
os.system('sed -i s/UNE1/'+str(potentialSteps[j])+'/g '+initial+'UNe')
os.system('sed -i s/ENDTIME/'+str(timeSteps[j*2])+'/g '+initial+'UNe')
os.system('sed -i s/UE/'+str(UNE_start*0.3)+'/g '+initial+'Ue')


#other initital conditions and BCs
os.system('sed -i s/ALPHAIN/'+str(ALPHAIN)+'/g '+initial+'alpha.gas')
os.system('sed -i s/CH2E/'+str(CH2E)+'/g '+initial+'C_H2.electrolyte')
os.system('sed -i s/CO2E/'+str(CO2E)+'/g '+initial+'C_O2.electrolyte')
os.system('sed -i s/CH2OE/'+str(CH2OE)+'/g '+initial+'C_H2O.electrolyte')
os.system('sed -i s/COHE/'+str(COHE)+'/g '+initial+'C_OH.electrolyte')
os.system('sed -i s/CH2G/'+str(CH2G)+'/g '+initial+'C_H2.gas')
os.system('sed -i s/CO2G/'+str(CO2G)+'/g '+initial+'C_O2.gas')
os.system('sed -i s/CH2OGAS/'+str(CO2G)+'/g '+initial+'C_H2O.gas')
os.system('sed -i s/EPSE/'+str(eps_E)+'/g '+initial+'eps')
os.system('sed -i s/EPSS/'+str(eps_S)+'/g '+initial+'eps')

#Case/constant/*Properties
constants=case+"/constant/transportProperties";
os.system('sed -i s/TAUC/'+str(TAUC)+'/g '+constants)
os.system('sed -i s/CH2SAT/'+str(CH2SAT)+'/g '+constants)
os.system('sed -i s/CO2SAT/'+str(CO2SAT)+'/g '+constants)
os.system('sed -i s/PNUM/'+str(PNUM)+'/g '+constants)
os.system('sed -i s/VapCase/'+VapourCase+'/g '+constants)
os.system('sed -i s/CMODEL/'+CModel+'/g '+constants)
os.system('sed -i s/NUMODEL/'+nuModel+'/g '+constants)
os.system('sed -i s/GG/'+str(GG)+'/g '+case+"/constant/g");
os.system('sed -i s/KPE/'+str(KPE)+'/g '+case+"/constant/fvOptions");
os.system('sed -i s/JLIM/'+str(JLim)+'/g '+constants);
os.system('sed -i s/ACPE/'+str(1-aPe)+'/g '+case+"/constant/reactionProperties");
os.system('sed -i s/AAPE/'+str(aPe)+'/g '+case+"/constant/reactionProperties");
os.system('sed -i s/AANE/'+str(1-aNe)+'/g '+case+"/constant/reactionProperties");
os.system('sed -i s/ACNE/'+str(aNe)+'/g '+case+"/constant/reactionProperties");
os.system('sed -i s/I0PE_REF/'+str(I0PE*I0_Scale)+'/g '+case+"/constant/reactionProperties");
os.system('sed -i s/I0NE_REF/'+str(I0NE*I0_Scale)+'/g '+case+"/constant/reactionProperties");
os.system('sed -i s/EAO_NE/'+str(Ea0_Ne)+'/g '+case+"/constant/reactionProperties");
os.system('sed -i s/EAO_PE/'+str(Ea0_Pe)+'/g '+case+"/constant/reactionProperties");
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#Use setFields to define fields.	
#prepare values to insert in setFields
#geometric values are defined at blockmesh
RPE	= r_p_E	 
RPS	= r_p_S	
EPSE = eps_E	
EPSS = eps_S	
	
setFieldsDictOrig=case+"/system/setFieldsDict.orig"
setFieldsDict=case+"/system/setFieldsDict"

os.system('sed -e s/RPE/'   + str(RPE)     + '/g  ' + setFieldsDictOrig+' > '+setFieldsDict)
#values
os.system('sed -i s/RPS/'      + str(RPS)  + '/g  ' + setFieldsDict)
os.system('sed -i s/EPSE/'      + str(EPSE)  + '/g  ' + setFieldsDict)
os.system('sed -i s/EPSS/'       + str(EPSS)   + '/g  ' + setFieldsDict)
os.system('sed -i s/TPE/'+str(TEL)+'/g  ' +setFieldsDict)
os.system('sed -i s/TNE/'+str(TEL)+'/g  ' +setFieldsDict)
os.system('sed -i s/TDI/'+str(TDI)+'/g  ' +setFieldsDict)
os.system('sed -i s/LC/'+str(LC)+'/g  ' +setFieldsDict)
os.system('sed -i s/WC/'+str(WC)+'/g  ' +setFieldsDict)
os.system('sed -i s/TAUS/'+str(tau_s)+'/g  ' +setFieldsDict)


os.system('echo running setFields...')
os.system('setFields -case '+ case+' > '+case+'/log.setFields')

	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	

	
	##########################################################################################
	
#Comment the following lines if you don't want to run the solver.        
decomposeParDictOrig=case+"/system/decomposeParDict.orig"
decomposeParDict=case+"/system/decomposeParDict"
os.system('sed -e s/NPROCESSOR/'+str(NProcessor)+'/g  '+decomposeParDictOrig+' > '+decomposeParDict)

##########################################################################################
if (int(NProcessor) > 1):
 os.system('echo '+str(NProcessor))
 decomposeParDictOrig=case+"/system/decomposeParDict.orig"
 decomposeParDict=case+"/system/decomposeParDict"
 os.system('sed -e s/NPROCESSOR/'+str(NProcessor)+'/g  '+decomposeParDictOrig+' > '+decomposeParDict)
 os.system('echo running decomposePar -force...')
 os.system('decomposePar -case '+ case+' -force > '+case+'/log.decomposePar')

 os.system('echo running: srun/mpiexec -np '+str(NProcessor)+' hisDriftFluxFoam -parallel -case '+ case+' > '+case+'/log.hisDriftFluxFoam ...')

 os.system('mpirun -np '+str(NProcessor)+' hisDriftFluxFoam -parallel -case '+ case+' > '+case+'/log.hisDriftFluxFoam')
 os.system('reconstructPar -case '+ case+' > '+case+'/log.ReconstructPar')
else:
 os.system('echo running hisDriftFluxFoam on single core')
 os.system('hisDriftFluxFoam -case '+ case+' > '+case+'/log.hisDriftFluxFoam')

set1Orig=case+"/system/set1.orig"
set1=case+"/system/set1"


os.system('sed -e s/XEnd/'   + str(LC*0.9*1e-3)     + '/g  ' + set1Orig+' > '+set1)
os.system('sed -i s/XMid/'+str(LC*0.5*1e-3)+'/g  ' +set1)
os.system('sed -i s/XStart/'+str(LC*0.1*1e-3)+'/g  ' +set1)
os.system('sed -i s/YTOP/'+str((TEL*2+TDI)*1e-3)+'/g  ' +set1)
os.system('sed -i s/ZMid/'+str((WC*0.5)*1e-3)+'/g  ' +set1)
os.system('echo running postProcess set1...')
os.system('postProcess -func set1 -case '+ case +' > '+ case+'/log.postProcessSet1')
os.system("python3 Cal_RMSE.py "+str(case)+" "+expL)
