import numpy as np
import sys
import os

#reading data from case
d=sys.argv[1]

f=open(d+'/postProcessing/PWall_Pe/0/surfaceFieldValue.dat',"r")
j=0
#loop to save log data.
simIRaw=[];
simURaw=[];
simI=[];
simU=[];
for line in f:
    if j>4:
    	simIRaw.append(abs(float(line.split()[2]))*(1000/(100*100)))#A/m² to mA/cm²
    j=j+1
f.close()
f=open(d+'/postProcessing/PWall_Ne/0/surfaceFieldValue.dat',"r")
j=0
for line in f:
    if j>4:
    	simURaw.append(abs(float(line.split()[4])))
    j=j+1
f.close()
for i in range(len(simURaw)):
	if i==(len(simURaw)-1):
		simI.append(round(simIRaw[i],3))
		simU.append(simURaw[i])
	elif simURaw[i]!=simURaw[i+1]:
		simI.append(round(simIRaw[i],3))
		simU.append(simURaw[i])
#readingReference
base="~/experimental_results/";
expL=sys.argv[2];
a,b=np.loadtxt(os.path.expanduser(base+expL), unpack=True)
DU=list(b)
DI=list(a)
print(DU)
#caclulating RMSE
SES=0	#squre error sum.
n=0
for i in range(len(DU)):
	for j in range(len(simU)):
		if DU[i]==simU[j]:
			SES=SES+pow(simI[j]-DI[i],2)
			n=n+1
RMSE=pow(SES/n,0.5)
print(RMSE)
#writting error to caseInformationFile
f=open(d+"/caseInformation.txt","r")
check=1;
check2=1;
for line in f:
	if line=="RMSE:\n":
		check=0;
	if (line.split(":")[0]=="relativeTo" and line.split(":")[1]==(base+expL+str("\n"))):
		print("study done before")
		check2=0;
f.close()
f=open(d+"/caseInformation.txt","a")
if check:	
	f.write("RMSE:\n")
	f.write("Error:"+str(RMSE)+str("\n"))
	f.write("relativeTo:"+base+expL+str("\n"))
elif check2:
	f.write("NthCheck:\n")
	f.write("Error:"+str(RMSE)+str("\n"))
	f.write("relativeTo:"+base+expL+str("\n"))
f.close()






