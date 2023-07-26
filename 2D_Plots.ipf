#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


function run_all()
start()
load_files()
twodim()
end

Function start()
variable m


///////// Files to load here //////////

variable numberofexperiments=8						// Change here for the number of experiments

make/o/T/n=(numberofexperiments) filelist	



filelist[0]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample22:"
filelist[1]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample23:"
filelist[2]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample24:"
filelist[3]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample25:"
filelist[4]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample26:"
filelist[5]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample27:"
filelist[6]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample28:"
filelist[7]="Macintosh HD:Users:Mathew:Documents:Current analysis:20230718_lysate:Sample29:"




end


function load_files()
wave/T filelist
string file="FRET_and_size.csv"
print file
variable a

for(a=0;a<(dimsize(filelist,0));a+=1)
string name=num2str(a)
newdatafolder/s $name

string string_to_load=filelist[a]+file
print string_to_load
LoadWave/J/D/W/K=0/A string_to_load
matrices_fret_run()
setdatafolder root:
endfor

end




Function matrices_fret_run()
wave FRET,Size
duplicate/o size,sizew
Make/O/N=(20,30) Matrix = 0
Make/O/N=(Dimsize(FRET,0),8) FretH
variable l,k,m,count1,count2
//Bin the size and FRET data
for(l=0;l<=(Dimsize(FRET,0));l+=1)
	FretH[l][0] = FRET[l]/0.05				// FRET Efficiency
	FretH[l][1] = SizeW[l]/5                 		// Size
	FretH[l][2] = round(FretH[l][0])			// Rounded FRET Efficiency
	FretH[l][3] = round(FretH[l][1])			// Rounded Size
	FretH[l][4] = FretH[l][0] - FretH[l][2]		// Needed to convert rounded to int
	FretH[l][5] = FretH[l][1] - FretH[l][3]		// Needed to convert rounded to int
if(FretH[l][4] < 0)
	FretH[l][6] = FretH[l][2] -1				// FRET Efficiency Int
else
	FretH[l][6] = FretH[l][2]				// FRET Efficiency Int
endif
if(FretH[l][5] < 0)
	FretH[l][7] = FretH[l][3] -1				// Size Int
else
	FretH[l][7] = FretH[l][3]				// Size Int
endif
endfor

////The following is data not accounting for Pascal

//Delete any data in which the size doesn't fit into the matrix
for(k=0;k<=(Dimsize(FretH,0));k+=1)
		if(FretH[k][6] > 20)
		DeletePoints k,1, FretH
		count1 +=1
		endif
endfor
for(k=0;k<=(Dimsize(FretH,0));k+=1)
	if(FretH[k][7] > 30)
	DeletePoints k,1, FretH
	count2 +=1
	endif
endfor
Printf "%d data points deleted due to FRET being too high.\r",count1
Printf "%d data points deleted due to Size being too high.\r",count2
make/o/n=1 toobig=count2

//Add the data to the matrix.
variable x,y
for(m=0;m<=(Dimsize(FretH,0));m+=1)
if(FretH[m][7]<=29)
if(FretH[m][6]<=19)
y = FretH[m][6]
x = FretH[m][7] 
duplicate/o matrix, matrix10
matrix[y][x]=matrix10[y][x]+1 // Adds 1 each time data needs to be added
endif
endif
endfor

////Accounting for Pascal
Make/O/N=(20,150) pascalm = 0
Make/O/N=(Dimsize(FRET,0),8) Pascal_Fret
for(l=0;l<=(Dimsize(FRET,0));l+=1)
	Pascal_Fret[l][0] = FRET[l]/0.05				// Bin FRET Efficiency
	Pascal_Fret[l][1] = SizeW[l]					//sizes not binned. 
	Pascal_Fret[l][2] = round(Pascal_Fret[l][0])
	Pascal_Fret[l][3] = round(Pascal_Fret[l][1])
	Pascal_Fret[l][4] = Pascal_Fret[l][0] - Pascal_Fret[l][2]
	Pascal_Fret[l][5] = Pascal_Fret[l][1] - Pascal_Fret[l][3]
if(Pascal_Fret[l][4] < 0)
	Pascal_Fret[l][6] = Pascal_Fret[l][2] -1
else
	Pascal_Fret[l][6] = Pascal_Fret[l][2]		// FRET Efficiency
endif
if(Pascal_FRET[l][5] < 0)
	Pascal_Fret[l][7] = Pascal_Fret[l][3] -1
else
	Pascal_Fret[l][7] = Pascal_Fret[l][3]		// Sizes
endif
endfor

//Delete points that won't fit into the 150x20 Matrix

for(k=0;k<=(Dimsize(Pascal_Fret,0));k+=1)
		if(Pascal_Fret[k][6] > 20)
		DeletePoints k,1, Pascal_Fret
		endif
endfor

for(k=0;k<=(Dimsize(Pascal_Fret,0));k+=1)
	if(pascal_FRET[k][7] > 150)
	DeletePoints k,1, Pascal_Fret
	endif

endfor

//Put data into the 150x20 Matrix

for(m=0;m<=(Dimsize(Pascal_Fret,0));m+=1)
variable count4
if(pascal_FRET[m][7]<=149)
if(pascal_FRET[m][6]<=19)
y = Pascal_Fret[m][6]
x = Pascal_Fret[m][7] 
duplicate/o pascalm, pascalm10
pascalm[y][x]=pascalm10[y][x]+1
endif
endif
endfor


//Multiply size distributions by Pascal factor to account for "invisible" oligomers

for(m=0;m<=(Dimsize(pascalm,0));m+=1)

// This accounts for pascal- I have commented out for the time being, since our size estimation is quite poor!

//pascalm[m][2]=pascalm[m][2]*2
//pascalm[m][3]=pascalm[m][3]*1.33
//pascalm[m][4]=pascalm[m][4]*1.14
//pascalm[m][5]=pascalm[m][5]*1.067



endfor


////Now to modify the matrix not accounted for Pascal
duplicate/o matrix, Matrix_pascal
for(m=0;m<=(Dimsize(Matrix_pascal,0));m+=1)

matrix_pascal[m][0]=pascalm[m][0]+pascalm[m][1]+pascalm[m][2]+pascalm[m][3]+pascalm[m][4]
matrix_pascal[m][1]=pascalm[m][5]+pascalm[m][6]+pascalm[m][7]+pascalm[m][8]+pascalm[m][9]

endfor
duplicate/o pascalm FRET_Matrix_1Bin
killwaves pascalm
killwaves pascalm10
//killwaves Pascal_FRET
killwaves matrix10
killwaves FretH
killwaves matrix
duplicate matrix_pascal FRET_Matrix




duplicate/o FRET_Matrix_1Bin FRET_Matrix_1Bin_Mass
variable n=1
for(k=0;k<=(Dimsize(FRET_Matrix_1Bin,1));k+=1)
for(m=0;m<=(Dimsize(FRET_Matrix_1Bin,0));m+=1)

FRET_Matrix_1Bin_Mass[m][k] = FRET_Matrix_1Bin[m][k]*k

n+=1
endfor
endfor
make/o/n=20 FRETaxis
make/o/n=150 size1axis
variable r = 0.025
variable t = 0
for(m=0;m<=19;m+=1)

FRETaxis[m] = r
r+=0.05

endfor
for(m=0;m<=149;m+=1)

size1axis[m] = t
t+=1



endfor
duplicate FRET_Matrix_1Bin_Mass,ln_FRET_Matrix_1Bin_Mass
variable a,b
for(a=0;a<150;a+=1)
for(b=0;b<20;b+=1)
if(ln(FRET_Matrix_1Bin_Mass[a][b])>0)
ln_FRET_Matrix_1Bin_Mass[a][b]=ln(FRET_Matrix_1Bin_Mass[a][b])
else
ln_FRET_Matrix_1Bin_Mass[a][b]=0
endif
endfor
endfor

Display;AppendMatrixContour ln_FRET_Matrix_1Bin_Mass vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,50
ModifyContour ln_FRET_Matrix_1Bin_Mass ctabLines={0,*,Rainbow256,1}
ModifyContour ln_FRET_Matrix_1Bin_Mass ctabLines={0,12,Rainbow256,1},autoLevels={0,12,100}
ModifyContour ln_FRET_Matrix_1Bin_Mass labels=0
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "FRET Efficiency"
Label left "Size of Oligomer"
ModifyGraph width=226.772,height=226.772
ColorScale/C/N=text0/F=0/A=RC/E contour=ln_FRET_Matrix_1Bin_Mass

duplicate/o matrix_pascal matrix_pascal_mass

for(m=0;m<=(Dimsize(matrix_pascal,0));m+=1)

matrix_pascal_mass[m][0] = matrix_pascal_mass[m][0]
matrix_pascal_mass[m][1] = matrix_pascal_mass[m][1]*5
matrix_pascal_mass[m][2] = matrix_pascal_mass[m][2]*10
matrix_pascal_mass[m][3] = matrix_pascal_mass[m][3]*15
matrix_pascal_mass[m][4] = matrix_pascal_mass[m][4]*20
matrix_pascal_mass[m][5] = matrix_pascal_mass[m][5]*25
matrix_pascal_mass[m][6] = matrix_pascal_mass[m][6]*30
matrix_pascal_mass[m][7] = matrix_pascal_mass[m][7]*35
matrix_pascal_mass[m][8] = matrix_pascal_mass[m][8]*40
matrix_pascal_mass[m][9] = matrix_pascal_mass[m][9]*45
matrix_pascal_mass[m][10] = matrix_pascal_mass[m][10]*50
matrix_pascal_mass[m][11] = matrix_pascal_mass[m][11]*55
matrix_pascal_mass[m][12] = matrix_pascal_mass[m][12]*60
matrix_pascal_mass[m][13] = matrix_pascal_mass[m][13]*65
matrix_pascal_mass[m][14] = matrix_pascal_mass[m][14]*70
matrix_pascal_mass[m][15] = matrix_pascal_mass[m][15]*75
matrix_pascal_mass[m][16] = matrix_pascal_mass[m][16]*80
matrix_pascal_mass[m][17] = matrix_pascal_mass[m][17]*85
matrix_pascal_mass[m][18] = matrix_pascal_mass[m][18]*90
matrix_pascal_mass[m][19] = matrix_pascal_mass[m][19]*95
matrix_pascal_mass[m][20] = matrix_pascal_mass[m][20]*100
matrix_pascal_mass[m][21] = matrix_pascal_mass[m][21]*105
matrix_pascal_mass[m][22] = matrix_pascal_mass[m][22]*110
matrix_pascal_mass[m][23] = matrix_pascal_mass[m][23]*115
matrix_pascal_mass[m][24] = matrix_pascal_mass[m][24]*120
matrix_pascal_mass[m][25] = matrix_pascal_mass[m][25]*125
matrix_pascal_mass[m][26] = matrix_pascal_mass[m][26]*130
matrix_pascal_mass[m][27] = matrix_pascal_mass[m][27]*135
matrix_pascal_mass[m][28] = matrix_pascal_mass[m][28]*140
matrix_pascal_mass[m][29] = matrix_pascal_mass[m][29]*145

endfor
make/o/n=30 size5axis
variable u=0
for(m=0;m<=29;m+=1)
size5axis[m] = u
u+=5
endfor

killwaves matrix_pascal
duplicate/o matrix_pascal_mass FRET_Mass_Matrix
killwaves matrix_pascal_mass

//variable fact=A_bursts[301]
for(x=0;x<(dimsize(FRET_Matrix_1Bin,0));x+=1)
for(y=0;y<(dimsize(FRET_Matrix_1Bin,1));y+=1)
FRET_Matrix_1bin[x][y]=FRET_Matrix_1bin[x][y]
FRET_Matrix_1bin_mass[x][y]=FRET_Matrix_1bin_mass[x][y]
endfor
endfor
End 





function copy_plots()
wave  FRET_Matrix_1Bin,fretaxis,size1axis

duplicate/o FRET_Matrix_1Bin,seelnFRET_Matrix_1Bin
variable a,b
for(a=0;a<150;a+=1)
for(b=0;b<20;b+=1)
if(ln(FRET_Matrix_1Bin[a][b])>0)
seelnFRET_Matrix_1Bin[a][b]=ln(FRET_Matrix_1Bin[a][b])
else
seelnFRET_Matrix_1Bin[a][b]=0
endif
endfor
endfor

end



function twodim()
kill()
setdatafolder root:
wave filelist
variable n
for(n=0;n<(dimsize(filelist,0));n+=1)
setdatafolder root:

string y=num2str(n)
setdatafolder $y
copy_plots()
variable m
wave seelnFRET_Matrix_1Bin

make/o/n=20 FRETaxis
make/o/n=150 size1axis
variable r = 0.025
variable t = 0
for(m=0;m<=19;m+=1)

FRETaxis[m] = r
r+=0.05

endfor
for(m=0;m<=149;m+=1)

size1axis[m] = t
t+=1



endfor

// scale etc
Display;AppendMatrixContour seelnFRET_Matrix_1Bin vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,150
ModifyContour seelnFRET_Matrix_1Bin labels=0
ModifyContour seelnFRET_Matrix_1Bin autoLevels={*,*,50}
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,12,Rainbow256,1},autoLevels={0,12,50}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
ModifyGraph manTick(left)={0,20,0,0},manMinor(left)={0,0}
ModifyGraph fSize(bottom)=8
ModifyGraph fSize(left)=8
ModifyGraph swapXY=1
ModifyGraph width=120
ModifyGraph height=100
ModifyGraph noLabel=2
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E "\Z08"+y
TextBox/C/N=text1/X=8.00/Y=0.00
ModifyContour seelnFRET_Matrix_1Bin autoLevels={0,11,50}
SetAxis bottom 1,50
ModifyGraph btLen=5,manTick(bottom)={0,10,0,0},manMinor(bottom)={0,0}
ModifyGraph mirror=1
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,11,Rainbow256,1}
ModifyContour seelnFRET_Matrix_1Bin autoLevels={0,10,50}
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,10,Rainbow256,1}





endfor
Display;AppendMatrixContour seelnFRET_Matrix_1Bin vs {FRETaxis,size1axis}
SetAxis bottom 0,1
SetAxis left 1,150
ModifyContour seelnFRET_Matrix_1Bin labels=0
ModifyContour seelnFRET_Matrix_1Bin autoLevels={*,*,50}
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,12,Rainbow256,1},autoLevels={0,12,50}
ModifyGraph manTick(left)={0,5,0,0},manMinor(left)={0,0}
Label bottom "\Z08 FRET Efficiency"
Label left "\Z08 Size of Oligomer"
ColorScale/C/N=text4/F=0/A=RC/E frame=0.00
TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=0.00/E "Scale bar etc"
ModifyGraph manTick(left)={0,20,0,0},manMinor(left)={0,0}
ModifyGraph fSize(bottom)=8
ModifyGraph fSize(left)=8
ModifyGraph swapXY=1
ModifyContour seelnFRET_Matrix_1Bin autoLevels={0,11,50}
SetAxis bottom 1,50
ModifyGraph btLen=5,manTick(bottom)={0,10,0,0},manMinor(bottom)={0,0}
ModifyGraph mirror=1
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,11,Rainbow256,1}
ModifyContour seelnFRET_Matrix_1Bin autoLevels={0,10,50}
ModifyContour seelnFRET_Matrix_1Bin ctabLines={0,10,Rainbow256,1}

newLayout/N=Histograms_1
AppendLayoutObject/r=(75,95,280,230)/F=0 graph graph0
AppendLayoutObject/r=(75,229,280,364)/F=0 graph graph2
AppendLayoutObject/r=(75,363,280,498)/F=0 graph graph4
AppendLayoutObject/r=(75,497,280,632)/F=0 graph graph6
AppendLayoutObject/r=(75,631,280,766)/F=0 graph graph8
AppendLayoutObject/r=(279,95,484,230)/F=0 graph graph1
AppendLayoutObject/r=(279,229,484,364)/F=0 graph graph3
AppendLayoutObject/r=(279,363,484,498)/F=0 graph graph5
AppendLayoutObject/r=(279,497,484,632)/F=0 graph graph7
AppendLayoutObject/r=(279,631,484,766)/F=0 graph graph9

newLayout/N=Histograms_2
AppendLayoutObject/r=(75,95,280,230)/F=0 graph graph10
AppendLayoutObject/r=(75,229,280,364)/F=0 graph graph12
AppendLayoutObject/r=(75,363,280,498)/F=0 graph graph14
AppendLayoutObject/r=(75,497,280,632)/F=0 graph graph16
AppendLayoutObject/r=(75,631,280,766)/F=0 graph graph18
AppendLayoutObject/r=(279,95,484,230)/F=0 graph graph11
AppendLayoutObject/r=(279,229,484,364)/F=0 graph graph13
AppendLayoutObject/r=(279,363,484,498)/F=0 graph graph15
AppendLayoutObject/r=(279,497,484,632)/F=0 graph graph17
AppendLayoutObject/r=(279,631,484,766)/F=0 graph graph19


//AppendLayoutObject graph graph2
end

function kill()

	variable                      winMask;
 
	variable                      i,n;
	variable                      all=0x1000+0x40+0x10+0x4+0x2+0x1;
	string                        theWins;
 
	winMask = !winMask ? all : winMask;
 
	theWins = winList("*",";","WIN:"+num2iStr(winMask & all));
	for(i=0,n=itemsInList(theWins,";") ; i<n ; i+=1)
		doWindow/K $stringFromList(i,theWins,";");
	endfor;
end