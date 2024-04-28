from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from abaqus import mdb, session
import os
import numpy as np
from math import log10

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
mdb.models['Model-1'].setValues(absoluteZero=-273.16, stefanBoltzmann=5.667e-11)

#Input_Data
h=800       #[Beam Height as 'mm']
b=300       #[Beam Width as 'mm']
L=3200      #[Beam Length as 'mm']

r1=10      #[Tie-Bar bar dia as 'mm']
r2=10     #[Top bar dia as 'mm']
r3=25      #[Bottom bar dia as 'mm']

bar_dia = [6, 8, 10, 12, 14, 16, 20, 22, 25, 28, 32, 40, 50]
bar_area = [28.3, 50.3, 78.5, 113, 154, 201, 314, 380, 491, 616, 804, 1257, 1964]

find_area = lambda dia: next((area for d, area in zip(bar_dia, bar_area) if d == dia), None)

A1=find_area(r1)        #[Tie-Bar bar area as 'mm2']
A2=find_area(r2)        #[Top Bar bar area as 'mm2']
A3=find_area(r3)        #[Bottom Bar bar area as 'mm2']

Cb=75    #Bottom Tie-Bar rebar face to beam surface distance 'mm'
Ct=38        #Top  Tie-Bar rebar face to beam surface distance 'mm'
Cs=Cb       #Side  Tie-Bar rebar face to beam surface distance 'mm'

LS=1000      #Loading positon distance
S1=100      #Support from origin
S2=L-S1    #Support apart from origin
ST=100     #Tie-Bar bar start
SP=150     #Tie-Bar bar spacing as 'mm'
N=(((L-2*ST)/SP)+1)      #NuMain-Barer of  Tie-Bar bar

L1=0.5*(L-LS)    #Load position 01 from origin
L2=0.5*(L+LS)      #Load position 02 from origin

FZ1=200         #Fire position 01 from origin
FZ2=L-FZ1         #Fire position 02 from origin

t=150       #Slab Thickness that cover the beam as 'mm'
S=h-t

#Time in minutes
t1=60      
t2=90
t3=120      
t4=150
t5=180
t6=0

m=25        #Mesh size as 'mm'

start_index = int(0.5 * L / m)
end_index = start_index + 1

Te= [20, 100, 200, 300, 400, 500, 600, 700, 800]    #Temperature Degree unit (Steel)
f = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9333, 0.8667, 0.8]   #Residual yield factor (Steel)
E = 200     #Steel Elastic Modulus unit GPa
Fy = 420   #Steel Yield Strength unit MPa
eu = 0.03       #Ultimate strain of Steel
Ve=0.3     #Poissons ratio steel

Tc = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]    #Temperature Degree unit (Concrete)
fc1 = [1.0, 1.0, 1.0, 1.0, 0.8486, 0.6956, 0.5426, 0.3896, 0.2366, 0.13, 0.05, 0.012]       #Concrete Compressive Factor  
fc2 = [1.0, 0.896, 0.766, 0.636, 0.506, 0.376, 0.246, 0.077, 0.062, 0.046, 0.031, 0.015]       #Concrete Tensile Factor

Fc = 30   #Concrete Compressive Strength unit MPa
fe = 0.33   #Elastic Factor
epc = [0.0025, 0.004, 0.0055, 0.007, 0.01, 0.015, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025] #peak strain
euc = [0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, 0.0475] #ultimate strain
Vec=0.2  #poissons ratio concrete

name = 'B-'+str(h)+'x'+str(b)+'-C'+str(Cb)

# For loop for Fire Curve
def calculate_values(T):
    results = ()
    for x in range(0, (T*60)+1, 60):
        # Calculate y for the current value of x
        y = 20 + (345 * log10(((8*x )/60)+ 1))
        # Create a tuple (x, y) and add it to the results tuple
        results += ((round(x, 5), round(y, 5)),)
    return results

#Material_Property_Create
# For loop for Steel Material
ms = ()
mse=()
for temp, rate in zip(Te, f):
    e = np.linspace(((Fy * rate) / (E * 1000)), eu, 10)
    St = Fy * rate * (1 + e)
    et = np.log(1 + e) - np.log((1+((Fy * rate) / (E * 1000))))  # Applying natural logarithm element-wise
    
    for s, e_val in zip(St, et):
        ms+=((round(s, 5), round(e_val, 5), temp),)

for temp in Te:
    mse += ((E*1000, Ve, temp),)

mdb.models['Model-1'].Material(name='Steel-Material')
mdb.models['Model-1'].materials['Steel-Material'].Density(table=((7.85e-09, ), ))
mdb.models['Model-1'].materials['Steel-Material'].Conductivity(table=((53.3, 
    20.0), (50.7, 100.0), (47.3, 200.0), (44.0, 300.0), (40.7, 400.0), (37.4, 
    500.0), (34.0, 600.0), (30.7, 700.0), (27.3, 800.0), (27.3, 900.0), (27.3, 
    1000.0), (27.3, 1100.0), (27.3, 1200.0)), temperatureDependency=ON)
mdb.models['Model-1'].materials['Steel-Material'].SpecificHeat(table=((
    439801760.0, 20.0), (487620000.0, 100.0), (529760000.0, 200.0), (
    564740000.0, 300.0), (605880000.0, 400.0), (666500000.0, 500.0), (
    760217391.3, 600.0), (1008157895.0, 700.0), (5000000000.0, 735.0), (
    803260869.6, 800.0), (650000000.0, 900.0), (650000000.0, 1000.0), (
    650000000.0, 1100.0), (650000000.0, 1200.0)), temperatureDependency=ON)
mdb.models['Model-1'].materials['Steel-Material'].Plastic(table=ms, temperatureDependency=ON)
mdb.models['Model-1'].materials['Steel-Material'].Elastic(table=mse, temperatureDependency=ON)

# For loop for Concrete-Material
mc = ()
mce = ()
mct = ()

# Loop through values of T and f simultaneously
for temp, rate, peak, ultimate in zip(Tc, fc1, epc, euc):
    # Calculate e for the current value of rate
    e = np.concatenate((np.linspace((0.5*rate*peak*fe), peak, 6, endpoint=False),
                        np.linspace(peak, ultimate, 10)))
    St = np.concatenate(([rate*Fc*fe], rate * Fc * ((2*(e[1:]/peak))/(1+(e[1:]/peak)**2))))
    et = np.concatenate(([0], e[1:]-((St[1:]*peak)/(2*Fc*rate)))) 
    # Flatten the nested tuples
    for s, e_val in zip(St, et):
        mc+=((round(s, 5), round(e_val, 5), temp),)

for temp, rate, peak in zip(Tc, fc1, epc):
    mce += ((round(((2*Fc*rate)/(peak)), 3), Vec, temp),)

for temp, rate,  rate1, peak in zip(Tc, fc2, fc1, epc):
    # Calculate e for the current value of rate
    e = np.array([((0.1*Fc*rate*peak)/(2*Fc*rate1)), 1.25*((0.1*Fc*rate*peak)/(2*Fc*rate1)),
                  4.0*((0.1*Fc*rate*peak)/(2*Fc*rate1)), 8.7*((0.1*Fc*rate*peak)/(2*Fc*rate1))])
    St = np.array([0.1*Fc*rate, 0.1*0.77*Fc*rate, 0.1*0.45*Fc*rate, 0.1*0.1*Fc*rate])
    et = e-((St*peak)/(2*Fc*rate1))  # Applying natural logarithm element-wise
    # Flatten the nested tuples
    for s, e_val in zip(St, et):
        mct+=((round(s, 5), round(e_val, 5), temp),)

mdb.models['Model-1'].Material(name='Concrete-Material')
mdb.models['Model-1'].materials['Concrete-Material'].Density(table=((2.4e-09, ), ))
mdb.models['Model-1'].materials['Concrete-Material'].Conductivity(table=((1.95, 
    20.0), (1.77, 100.0), (1.55, 200.0), (1.36, 300.0), (1.19, 400.0), (1.04, 
    500.0), (0.91, 600.0), (0.81, 700.0), (0.72, 800.0), (0.66, 900.0), (0.62, 
    1000.0), (0.6, 1100.0), (0.6, 1200.0)), temperatureDependency=ON)
mdb.models['Model-1'].materials['Concrete-Material'].Elastic(table=mce, temperatureDependency=ON)
mdb.models['Model-1'].materials['Concrete-Material'].SpecificHeat(table=((
    900000000.0, 20.0), (900000000.0, 100.0), (1000000000.0, 200.0), (
    1050000000.0, 300.0), (1100000000.0, 400.0), (1100000000.0, 500.0), (
    1100000000.0, 600.0), (1100000000.0, 700.0), (1100000000.0, 800.0), (
    1100000000.0, 900.0), (1100000000.0, 1000.0), (1100000000.0, 1100.0), (
    1100000000.0, 1200.0)), temperatureDependency=ON)
mdb.models['Model-1'].materials['Concrete-Material'].ConcreteDamagedPlasticity(
    table=((31.0, 0.1, 1.16, 0.667, 0.0001), ))
mdb.models['Model-1'].materials['Concrete-Material'].concreteDamagedPlasticity.ConcreteCompressionHardening(
    table=mc, temperatureDependency=ON)
mdb.models['Model-1'].materials['Concrete-Material'].concreteDamagedPlasticity.ConcreteTensionStiffening(
    table=mct, temperatureDependency=ON)

#All types of data input end here
'''***************************'''

# Change these variables to your desired folder path
folder_path = r'F:\MSc\Post-Fire-'+str(name)
model_name = str(name)

# Check if the folder exists, if not, create it
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# Change directory to the specified folder
os.chdir(folder_path)

# Save the model in the specified folder
mdb.saveAs(pathName=os.path.join(folder_path, model_name + '.cae'))

#Create Part by COORDINATE (Only Main Body Created)
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0, 0), point2=(b, h))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Beam', type=DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Beam'].BaseSolidExtrude(depth=L, sketch=mdb.models['Model-1'].sketches['__profile__'])

#Create Datum Plane
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=S1, principalPlane=XYPLANE)      #This line mean datum[2]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=FZ1, principalPlane=XYPLANE)      #This line mean datum[3]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=L1, principalPlane=XYPLANE)      #This line mean datum[4]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=0.5*L, principalPlane=XYPLANE)      #This line mean datum[5]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=L2, principalPlane=XYPLANE)      #This line mean datum[6]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=FZ2, principalPlane=XYPLANE)      #This line mean datum[7]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=S2, principalPlane=XYPLANE)      #This line mean datum[8]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=0.5*b, principalPlane=YZPLANE)      #This line mean datum[9]
mdb.models['Model-1'].parts['Beam'].DatumPlaneByPrincipalPlane(offset=S, principalPlane=XZPLANE)      #This line mean datum[10]

#Create Partition
#Use loop to avoid multiple line code
for i in range(2, 9):
        mdb.models['Model-1'].parts['Beam'].PartitionCellByDatumPlane(cells=
        mdb.models['Model-1'].parts['Beam'].cells.findAt(((0.5*b, h, S2+5), )), datumPlane=
        mdb.models['Model-1'].parts['Beam'].datums[i])

#Create Partition
#From Here No loop is used
mdb.models['Model-1'].parts['Beam'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Beam'].cells.findAt(((b, 0, S1-5), ), ((b, 0, FZ1-5), ), ((b, 0, L1-5), ), ((b, 0, L1+5), ), ((b, 0, L2-5), ), ((b, 0, FZ2-5), ), ((b, 0, S2-5), ), ((b, 0, S2+5), ), ), 
    datumPlane=mdb.models['Model-1'].parts['Beam'].datums[9])

mdb.models['Model-1'].parts['Beam'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Beam'].cells.findAt(((b, 0, S1-5), ), ((b, 0, FZ1-5), ), ((b, 0, L1-5), ), ((b, 0, L1+5), ), ((b, 0, L2-5), ), ((b, 0, FZ2-5), ), ((b, 0, S2-5), ), ((b, 0, S2+5), ), 
    ((0, 0, S1-5), ), ((0, 0, FZ1-5), ), ((0, 0, L1-5), ), ((0, 0, L1+5), ), ((0, 0, L2-5), ), ((0, 0, FZ2-5), ), ((0, 0, S2-5), ), ((0, 0, S2+5), ),  ), 
    datumPlane=mdb.models['Model-1'].parts['Beam'].datums[10])

#Create Part by COORDINATE (Remain Parts Created)
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(L, 0.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.5*L, 0.0))
mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(addUndoState=False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.5*L, 0.0), ))

mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Top-Bar', type=DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Top-Bar'].BaseWire(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].Part(name='Main-Bar', objectToCopy=mdb.models['Model-1'].parts['Top-Bar'])

mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2000.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0, 0), point2=(b-2*Cs-r1, h-Ct-Cb-r1))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Tie-Bar', type=DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Tie-Bar'].BaseWire(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

#Assembly All Parts

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Top-Bar-1', part=mdb.models['Model-1'].parts['Top-Bar'])
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 0.0), direction2=(0.0, 1.0, 0.0), instanceList=('Top-Bar-1', ),
 number1=1, number2=2, spacing1=L, spacing2=h-Ct-Cb-2*r1-0.5*r2-0.5*r3)
del mdb.models['Model-1'].rootAssembly.features['Top-Bar-1']

mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Main-Bar-1', part=mdb.models['Model-1'].parts['Main-Bar'])
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 0.0), direction2=(0.0, 0.0, 1.0), instanceList=('Top-Bar-1-lin-1-2', 'Main-Bar-1'),
 number1=1, number2=2, spacing1=L, spacing2=b-2*Cs-2*r1-r3)
 
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 10.0, 0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Top-Bar-1-lin-1-2', 'Main-Bar-1', 'Top-Bar-1-lin-1-2-lin-1-2', 'Main-Bar-1-lin-1-2'))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Tie-Bar-1', part=mdb.models['Model-1'].parts['Tie-Bar'])

mdb.models['Model-1'].rootAssembly.translate(instanceList=('Top-Bar-1-lin-1-2', ), vector=(-0.5*(r3-r2), 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Top-Bar-1-lin-1-2-lin-1-2', ), vector=(0.5*(r3-r2), 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Tie-Bar-1', ), vector=(-0.5*(r3+r1), -0.5*(r3+r1), 0.0))

mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 0.0), direction2=(0.0, 0.0, -1.0), instanceList=('Tie-Bar-1', ), number1=1, number2=2, spacing1=112.0, spacing2=ST)

mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 0.0), direction2=(0.0, 0.0, -1.0), instanceList=('Tie-Bar-1-lin-1-2', ), 
    number1=1, number2=N, spacing1=112.0, spacing2=SP)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Beam-1', part=mdb.models['Model-1'].parts['Beam'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Beam-1', ), vector=(-(0.5*r3+r1+Cs), -(0.5*r3+r1+Cb), -L))
del mdb.models['Model-1'].rootAssembly.features['Tie-Bar-1']

mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 1.0, 0.0), instanceList=('Main-Bar-1', ), number1=2, 
    number2=1, spacing1=0.5*(b-2*Cs-2*r1-r3), spacing2=1.0)

#Material_Property_Related_Step
 
mdb.models['Model-1'].HomogeneousSolidSection(material='Concrete-Material', name='Beam', thickness=None)
mdb.models['Model-1'].TrussSection(area=A3, material='Steel-Material', name='Main-Bar')
mdb.models['Model-1'].TrussSection(area=A2, material='Steel-Material', name='Top-Bar')
mdb.models['Model-1'].TrussSection(area=A1, material='Steel-Material', name='Tie-Bar')

#Beam-Set-Create
mdb.models['Model-1'].parts['Beam'].Set(cells=
    mdb.models['Model-1'].parts['Beam'].cells.findAt(
    ((0.5*b-5, h, S1-5), ), ((0.5*b+5, h, S1-5), ), ((0.5*b-5, 0, S1-5), ), ((0.5*b+5, 0, S1-5), ), 
    ((0.5*b-5, h, S1+5), ), ((0.5*b+5, h, S1+5), ), ((0.5*b-5, 0, S1+5), ), ((0.5*b+5, 0, S1+5), ),
    ((0.5*b-5, h, S2-5), ), ((0.5*b+5, h, S2-5), ), ((0.5*b-5, 0, S2-5), ), ((0.5*b+5, 0, S2-5), ), 
    ((0.5*b-5, h, S2+5), ), ((0.5*b+5, h, S2+5), ), ((0.5*b-5, 0, S2+5), ), ((0.5*b+5, 0, S2+5), ),
    ((0.5*b-5, h, L1-5), ), ((0.5*b+5, h, L1-5), ), ((0.5*b-5, 0, L1-5), ), ((0.5*b+5, 0, L1-5), ), 
    ((0.5*b-5, h, L1+5), ), ((0.5*b+5, h, L1+5), ), ((0.5*b-5, 0, L1+5), ), ((0.5*b+5, 0, L1+5), ),
    ((0.5*b-5, h, L2-5), ), ((0.5*b+5, h, L2-5), ), ((0.5*b-5, 0, L2-5), ), ((0.5*b+5, 0, L2-5), ), 
    ((0.5*b-5, h, L2+5), ), ((0.5*b+5, h, L2+5), ), ((0.5*b-5, 0, L2+5), ), ((0.5*b+5, 0, L2+5), ),
    ), name='Set-Beam')

#Material-Assign
mdb.models['Model-1'].parts['Beam'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Beam'].sets['Set-Beam'], sectionName='Beam', thicknessAssignment=FROM_SECTION)
    
mdb.models['Model-1'].parts['Main-Bar'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    edges=mdb.models['Model-1'].parts['Main-Bar'].edges.findAt(((0.5*L, 0.0, 0.0), ), )), sectionName='Main-Bar', thicknessAssignment=FROM_SECTION)

mdb.models['Model-1'].parts['Top-Bar'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    edges=mdb.models['Model-1'].parts['Top-Bar'].edges.findAt(((0.5*L, 0.0, 0.0), ), )), sectionName='Top-Bar', thicknessAssignment=FROM_SECTION)

mdb.models['Model-1'].parts['Tie-Bar'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    edges=mdb.models['Model-1'].parts['Tie-Bar'].edges.findAt(((0.5*(b-2*Cs-r1), 0.0, 0.0), ), ((0.5*(b-2*Cs-r1), h-Ct-Cb-r1, 0.0), ), 
    ((0.0, 0.5*(h-Ct-Cb-r1), 0.0), ), ((b-2*Cs-r1, 0.5*(h-Ct-Cb-r1), 0.0), ), )), sectionName='Tie-Bar', thicknessAssignment=FROM_SECTION)

#Step_&_Field_Output_Create
mdb.models['Model-1'].rootAssembly.Set(name='U', vertices=
    mdb.models['Model-1'].rootAssembly.instances['Beam-1'].vertices.findAt(((0.5*b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb), -0.5*L), )))
    
mdb.models['Model-1'].CoupledTempDisplacementStep(deltmx=25.0, initialInc=60.0, maxInc=300.0, maxNumInc=1000000, minInc=0.05, name='Fire', previous='Initial', timePeriod=(t1*60))
mdb.models['Model-1'].StaticStep(initialInc=0.1, maxInc=0.2, maxNumInc=1000000, minInc=1e-15, name='PF-Load', previous='Fire')
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=('NT', ))
mdb.models['Model-1'].FieldOutputRequest(createStepName='PF-Load', name='F-Output-2', variables=('PE', 'VEEQ', 'U', 'RF'))

#Predefined-Field-Temperature
mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(20.0, ), name='Predefined Field-1', region=
    mdb.models['Model-1'].rootAssembly.instances['Beam-1'].sets['Set-Beam'])

#Fire-Surface-Create
mdb.models['Model-1'].rootAssembly.Surface(name='Fire', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Beam-1'].faces.findAt(
    ((b-(0.5*r3+r1+Cs), 0, -FZ2+5), ), ((-(0.5*r3+r1+Cs), 0, -FZ2+5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb), -FZ2+5), ),  ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb), -FZ2+5), ), 
    ((b-(0.5*r3+r1+Cs), 0, -FZ1-5), ), ((-(0.5*r3+r1+Cs), 0, -FZ1-5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb), -FZ1-5), ),  ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb), -FZ1-5), ), 
    ((b-(0.5*r3+r1+Cs), 0, -0.5*L+5), ), ((-(0.5*r3+r1+Cs), 0, -0.5*L+5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb), -0.5*L+5), ),  ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb), -0.5*L+5), ), 
    ((b-(0.5*r3+r1+Cs), 0, -0.5*L-5), ), ((-(0.5*r3+r1+Cs), 0, -0.5*L-5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb), -0.5*L-5), ),  ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb), -0.5*L-5), ), 
      ))
 
mdb.models['Model-1'].rootAssembly.Surface(name='Slab', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Beam-1'].faces.findAt(
    ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -FZ2+5), ), ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -FZ2+5), ), 
    ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -FZ1-5), ), ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -FZ1-5), ), 
    ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -0.5*L+5), ), ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -0.5*L+5), ), 
    ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -0.5*L-5), ), ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -0.5*L-5), ), 
      ))
 
mdb.models['Model-1'].rootAssembly.Surface(name='NoFire', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Beam-1'].faces.findAt(
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -FZ2+5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -FZ2+5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -FZ1-5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -FZ1-5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -0.5*L+5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -0.5*L+5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -0.5*L-5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -0.5*L-5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -S1+5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -S2+5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -S1-5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -S2-5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -S1+5), ), ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -S2+5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb), -S1-5), ), ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb), -S2-5), ), 
    ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -S1+5), ), ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -S2+5), ), 
    ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -S1-5), ), ((b-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5, -S2-5), ), 
    ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5,  -S1+5), ), ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5,  -S2+5), ), 
    ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5,  -S1-5), ), ((-(0.5*r3+r1+Cs), h-(0.5*r3+r1+Cb)-5,  -S2-5), ), 
    ((-(0.5*r3+r1+Cs), -(0.5*r3+r1+Cb)+5,  -S1+5), ), ((-(0.5*r3+r1+Cs), -(0.5*r3+r1+Cb)+5,  -S2+5), ), 
    ((-(0.5*r3+r1+Cs), -(0.5*r3+r1+Cb)+5,  -S1-5), ), ((-(0.5*r3+r1+Cs), -(0.5*r3+r1+Cb)+5,  -S2-5), ), 
	((b-(0.5*r3+r1+Cs), 0,  -S1+5), ), ((b-(0.5*r3+r1+Cs), 0,  -S2+5), ), 
    ((b-(0.5*r3+r1+Cs), 0,  -S1-5), ), ((b-(0.5*r3+r1+Cs), 0, -S2-5), ), 
	((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb),  -S1+5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb),  -S2+5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb),  -S1-5), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb), -S2-5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb),  -S1+5), ), ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb),  -S2+5), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb),  -S1-5), ), ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb), -S2-5), ),
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb)-5,  0), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb)-5,  0), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb)+5,  0), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb)+5,  0), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, h-(0.5*r3+r1+Cb)-5,  -L), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-(0.5*r3+r1+Cb)-5,  -L), ), 
    ((0.5*b-(0.5*r3+r1+Cs)-5, -(0.5*r3+r1+Cb)+5,  -L), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -(0.5*r3+r1+Cb)+5,  -L), ),       
    ))
    
 #Fire-Property-Impliment
mdb.models['Model-1'].FilmConditionProp(dependencies=0, name='NF', property=((
    0.009, ), ), temperatureDependency=OFF)
mdb.models['Model-1'].FilmConditionProp(dependencies=0, name='F', property=((
    0.025, ), ), temperatureDependency=OFF)
mdb.models['Model-1'].TabularAmplitude(data=calculate_values(t1), name='F-'+str(t1), smooth=SOLVER_DEFAULT, 
    timeSpan=STEP)
mdb.models['Model-1'].FilmCondition(createStepName='Fire', definition=
    PROPERTY_REF, interactionProperty='F', name='Int-1', sinkAmplitude='F-'+str(t1), 
    sinkDistributionType=UNIFORM, sinkFieldName='', sinkTemperature=1.0, 
    surface=mdb.models['Model-1'].rootAssembly.surfaces['Fire'])
mdb.models['Model-1'].interactions.changeKey(fromName='Int-1', toName='F')
mdb.models['Model-1'].FilmCondition(createStepName='Fire', definition=
    PROPERTY_REF, interactionProperty='F', name='Slab', sinkAmplitude='F-'+str(t1), 
    sinkDistributionType=UNIFORM, sinkFieldName='', sinkTemperature=1.0, 
    surface=mdb.models['Model-1'].rootAssembly.surfaces['Slab'])
mdb.models['Model-1'].FilmCondition(createStepName='Fire', definition=
    PROPERTY_REF, interactionProperty='NF', name='NF', sinkAmplitude='F-'+str(t1), 
    sinkDistributionType=UNIFORM, sinkFieldName='', sinkTemperature=1.0, 
    surface=mdb.models['Model-1'].rootAssembly.surfaces['NoFire'])
mdb.models['Model-1'].RadiationToAmbient(ambientTemperature=1.0, 
    ambientTemperatureAmp='F-'+str(t1), createStepName='Fire', distributionType=
    UNIFORM, emissivity=0.8, field='', name='Fire', radiationType=AMBIENT, 
    surface=mdb.models['Model-1'].rootAssembly.surfaces['Fire'])

#Embedded-Region

mdb.models['Model-1'].rootAssembly.regenerate()

edges_list = (
    mdb.models['Model-1'].rootAssembly.instances['Top-Bar-1-lin-1-2'].edges.getSequenceFromMask(mask=('[#1]', ), ) +
    mdb.models['Model-1'].rootAssembly.instances['Main-Bar-1'].edges.getSequenceFromMask(mask=('[#1]', ), ) +
    mdb.models['Model-1'].rootAssembly.instances['Top-Bar-1-lin-1-2-lin-1-2'].edges.getSequenceFromMask(mask=('[#1]', ), ) +
    mdb.models['Model-1'].rootAssembly.instances['Main-Bar-1-lin-1-2'].edges.getSequenceFromMask(mask=('[#1]', ), ) +
    mdb.models['Model-1'].rootAssembly.instances['Tie-Bar-1-lin-1-2'].edges.getSequenceFromMask(mask=('[#f]', ), )
)

for i in range(N-1):
    edges_list += mdb.models['Model-1'].rootAssembly.instances['Tie-Bar-1-lin-1-2-lin-1-' + str(i+2)].edges.getSequenceFromMask(mask=('[#f]', ), )

edges_list += mdb.models['Model-1'].rootAssembly.instances['Main-Bar-1-lin-2-1'].edges.getSequenceFromMask(mask=('[#1]', ), )

mdb.models['Model-1'].EmbeddedRegion(
    absoluteTolerance=0.0, 
    embeddedRegion=Region(edges=edges_list),
    fractionalTolerance=0.05, 
    hostRegion=None, 
    name='Embedded-All', 
    toleranceMethod=BOTH, 
    weightFactorTolerance=1e-06
)

#Reference_Point_Create
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.5*(b-r3)-r1-Cs, h-0.5*r3-r1-Cb, -L1))
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.5*(b-r3)-r1-Cs, h-0.5*r3-r1-Cb, -L2))
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.5*(b-r3)-r1-Cs, -0.5*r3-r1-Cb, -S1))
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.5*(b-r3)-r1-Cs, -0.5*r3-r1-Cb, -S2))

mdb.models['Model-1'].rootAssembly.Set(referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints.findAt((0.5*(b-r3)-r1-Cs, h-0.5*r3-r1-Cb, -L1),), ), name='RP1')
mdb.models['Model-1'].rootAssembly.Set(referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints.findAt((0.5*(b-r3)-r1-Cs, h-0.5*r3-r1-Cb, -L2),), ), name='RP2')
mdb.models['Model-1'].rootAssembly.Set(referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints.findAt((0.5*(b-r3)-r1-Cs, -0.5*r3-r1-Cb, -S1),), ), name='RP3')
mdb.models['Model-1'].rootAssembly.Set(referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints.findAt((0.5*(b-r3)-r1-Cs, -0.5*r3-r1-Cb, -S2),), ), name='RP4')

#Rigid_Body_Tie_Create
mdb.models['Model-1'].RigidBody(name='L1', refPointRegion=
    mdb.models['Model-1'].rootAssembly.sets['RP1'], tieRegion=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['Beam-1'].edges.findAt((
    (0.5*b-(0.5*r3+r1+Cs)-5, h-0.5*r3-r1-Cb, -L1), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-0.5*r3-r1-Cb, -L1), ), )))
mdb.models['Model-1'].RigidBody(name='L2', refPointRegion=
    mdb.models['Model-1'].rootAssembly.sets['RP2'], tieRegion=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['Beam-1'].edges.findAt((
    (0.5*b-(0.5*r3+r1+Cs)-5, h-0.5*r3-r1-Cb, -L2), ), ((0.5*b-(0.5*r3+r1+Cs)+5, h-0.5*r3-r1-Cb, -L2), ), )))
mdb.models['Model-1'].RigidBody(name='S1', refPointRegion=
    mdb.models['Model-1'].rootAssembly.sets['RP3'], tieRegion=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['Beam-1'].edges.findAt((
    (0.5*b-(0.5*r3+r1+Cs)-5, -0.5*r3-r1-Cb, -S1), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -0.5*r3-r1-Cb, -S1), ), )))
mdb.models['Model-1'].RigidBody(name='S2', refPointRegion=
    mdb.models['Model-1'].rootAssembly.sets['RP4'], tieRegion=Region(
    edges=mdb.models['Model-1'].rootAssembly.instances['Beam-1'].edges.findAt((
    (0.5*b-(0.5*r3+r1+Cs)-5, -0.5*r3-r1-Cb, -S2), ), ((0.5*b-(0.5*r3+r1+Cs)+5, -0.5*r3-r1-Cb, -S2), ), )))

#Boundary_Condition_Create
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-L', 
    region=mdb.models['Model-1'].rootAssembly.sets['RP3'], u1=SET, u2=SET, u3=
    SET, ur1=UNSET, ur2=SET, ur3=SET)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-R', 
    region=mdb.models['Model-1'].rootAssembly.sets['RP4'], u1=SET, u2=SET, 
    u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='PF-Load', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'UL1', region=mdb.models['Model-1'].rootAssembly.sets['RP1'], u1=UNSET, u2=
    -0.2*h, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='PF-Load', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'UL2', region=mdb.models['Model-1'].rootAssembly.sets['RP2'], u1=UNSET, u2=
    -0.2*h, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)


#Mesh_Create
mdb.models['Model-1'].parts['Beam'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=m)
mdb.models['Model-1'].parts['Beam'].generateMesh()
mdb.models['Model-1'].parts['Beam'].setElementType(elemTypes=(ElemType(
    elemCode=C3D8T, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    distortionControl=DEFAULT), ElemType(elemCode=C3D6T, elemLibrary=STANDARD), 
    ElemType(elemCode=C3D4T, elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['Beam'].cells.findAt(
    ((0.5*b-5, h, S1-5), ), ((0.5*b+5, h, S1-5), ), ((0.5*b-5, 0, S1-5), ), ((0.5*b+5, 0, S1-5), ), 
    ((0.5*b-5, h, S1+5), ), ((0.5*b+5, h, S1+5), ), ((0.5*b-5, 0, S1+5), ), ((0.5*b+5, 0, S1+5), ),
    ((0.5*b-5, h, S2-5), ), ((0.5*b+5, h, S2-5), ), ((0.5*b-5, 0, S2-5), ), ((0.5*b+5, 0, S2-5), ), 
    ((0.5*b-5, h, S2+5), ), ((0.5*b+5, h, S2+5), ), ((0.5*b-5, 0, S2+5), ), ((0.5*b+5, 0, S2+5), ),
    ((0.5*b-5, h, L1-5), ), ((0.5*b+5, h, L1-5), ), ((0.5*b-5, 0, L1-5), ), ((0.5*b+5, 0, L1-5), ), 
    ((0.5*b-5, h, L1+5), ), ((0.5*b+5, h, L1+5), ), ((0.5*b-5, 0, L1+5), ), ((0.5*b+5, 0, L1+5), ),
    ((0.5*b-5, h, L2-5), ), ((0.5*b+5, h, L2-5), ), ((0.5*b-5, 0, L2-5), ), ((0.5*b+5, 0, L2-5), ), 
    ((0.5*b-5, h, L2+5), ), ((0.5*b+5, h, L2+5), ), ((0.5*b-5, 0, L2+5), ), ((0.5*b+5, 0, L2+5), ),
    ), ))

mdb.models['Model-1'].parts['Main-Bar'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=m)
mdb.models['Model-1'].parts['Main-Bar'].generateMesh()
mdb.models['Model-1'].parts['Main-Bar'].setElementType(elemTypes=(ElemType(
    elemCode=T3D2T, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].parts['Main-Bar'].edges.findAt(((0.5*L, 0.0, 0.0), )), ))
    
mdb.models['Model-1'].parts['Tie-Bar'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=m)
mdb.models['Model-1'].parts['Tie-Bar'].generateMesh()
mdb.models['Model-1'].parts['Tie-Bar'].setElementType(elemTypes=(ElemType(
    elemCode=T3D2T, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].parts['Tie-Bar'].edges.findAt(((0.5*(b-2*Cs-r1), 0.0, 0.0), ), ((0.5*(b-2*Cs-r1), h-Ct-Cb-r1, 0.0), ), 
    ((0.0, 0.5*(h-Ct-Cb-r1), 0.0), ), ((b-2*Cs-r1, 0.5*(h-Ct-Cb-r1), 0.0), ), ), ))
    
mdb.models['Model-1'].parts['Top-Bar'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=m)
mdb.models['Model-1'].parts['Top-Bar'].generateMesh()
mdb.models['Model-1'].parts['Top-Bar'].setElementType(elemTypes=(ElemType(
    elemCode=T3D2T, elemLibrary=STANDARD), ), regions=(
    mdb.models['Model-1'].parts['Top-Bar'].edges.findAt(((0.5*L, 0.0, 0.0), )), ))

#History_Output_Request
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.DatumPointByCoordinate(coords=(b-2*Cs-2*r1-r3, 0.0, -0.5*L))
mdb.models['Model-1'].rootAssembly.Set(name='TC-R', nodes=
    mdb.models['Model-1'].rootAssembly.instances['Main-Bar-1-lin-1-2'].nodes[start_index:end_index])
del mdb.models['Model-1'].historyOutputRequests['H-Output-1']
mdb.models['Model-1'].HistoryOutputRequest(createStepName='PF-Load', 
    name='Displacement', rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.sets['U'], sectionPoints=DEFAULT, 
    variables=('U2', ))
mdb.models['Model-1'].HistoryOutputRequest(createStepName='Fire', name=
    'Temperature', rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.sets['TC-R'], sectionPoints=
    DEFAULT, variables=('NT', ))
mdb.models['Model-1'].historyOutputRequests['Temperature'].deactivate(
    'PF-Load')
mdb.models['Model-1'].HistoryOutputRequest(createStepName='PF-Load', 
    name='Force', rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.sets['RP3'], sectionPoints=
    DEFAULT, variables=('RF2', ))

#Model_Name_Change
mdb.models.changeKey(fromName='Model-1', toName=str(name)+'-F-'+str(t1))

#Remain model create
mdb.Model(name=str(name)+'-F-'+str(t2), objectToCopy=mdb.models[str(name)+'-F-'+str(t1)])
mdb.Model(name=str(name)+'-F-'+str(t3), objectToCopy=mdb.models[str(name)+'-F-'+str(t1)])
mdb.Model(name=str(name)+'-F-'+str(t4), objectToCopy=mdb.models[str(name)+'-F-'+str(t1)])
mdb.Model(name=str(name)+'-F-'+str(t5), objectToCopy=mdb.models[str(name)+'-F-'+str(t1)])
mdb.Model(name=str(name)+'-F-'+str(t6), objectToCopy=mdb.models[str(name)+'-F-'+str(t1)])

t = [t2, t3, t4, t5]  # Assuming t2, t3, ..., t8 are defined elsewhere

for i in range(2, 6):
    mdb.models[str(name)+'-F-'+str(t[i-2])].TabularAmplitude(data=calculate_values(t[i-2]), name='F-'+str(t[i-2]), smooth=SOLVER_DEFAULT, timeSpan=STEP)
    mdb.models[str(name)+'-F-'+str(t[i-2])].steps['Fire'].setValues(timePeriod=(t[i-2]*60))
    mdb.models[str(name)+'-F-'+str(t[i-2])].interactions['F'].setValues(definition=PROPERTY_REF
        , interactionProperty='F', sinkAmplitude='F-'+str(t[i-2]), sinkTemperature=1.0)
    mdb.models[str(name)+'-F-'+str(t[i-2])].interactions['Fire'].setValues(ambientTemperature=
        1.0, ambientTemperatureAmp='F-'+str(t[i-2]), distributionType=UNIFORM, emissivity=0.8, field='', radiationType=AMBIENT)
    mdb.models[str(name)+'-F-'+str(t[i-2])].interactions['NF'].setValues(definition=
        PROPERTY_REF, interactionProperty='NF', sinkAmplitude='F-'+str(t[i-2]), sinkTemperature=1.0)
    mdb.models[str(name)+'-F-'+str(t[i-2])].interactions['Slab'].setValues(definition=
        PROPERTY_REF, interactionProperty='F', sinkAmplitude='F-'+str(t[i-2]), sinkTemperature=1.0)

#Ambient Model Create
del mdb.models[str(name)+'-F-'+str(t6)].steps['Fire']

ms20 = ()
e = np.linspace(((Fy) / (E * 1000)), eu, 10)
St = Fy * (1 + e)
et = np.log(1 + e) - np.log((1+((Fy) / (E * 1000))))

for s, e_val in zip(St, et):
    ms20+=((round(s, 5), round(e_val, 5)),)

mdb.models[str(name)+'-F-'+str(t6)].Material(name='Steel-Material-20')
mdb.models[str(name)+'-F-'+str(t6)].materials['Steel-Material-20'].Density(table=((7.85e-09, ), ))
mdb.models[str(name)+'-F-'+str(t6)].materials['Steel-Material-20'].Plastic(table=ms20, temperatureDependency=OFF)
mdb.models[str(name)+'-F-'+str(t6)].materials['Steel-Material-20'].Elastic(table=((E*1000, Ve), ), temperatureDependency=OFF)

epc20 = 0.0025
euc20 = 0.02
mce20 = ()
mc20 = ()
mct20 = ()

e = np.concatenate((np.linspace((0.5*epc20*fe), epc20, 6, endpoint=False), np.linspace(epc20, euc20, 10)))
St = np.concatenate(([Fc*fe], Fc * ((2*(e[1:]/epc20))/(1+(e[1:]/epc20)**2))))
et = np.concatenate(([0], e[1:]-((St[1:]*epc20)/(2*Fc)))) 

for s, e_val in zip(St, et):
        mc20 +=((round(s, 5), round(e_val, 5)),)

mce20 += ((round(((2*Fc)/(epc20)), 3), Vec),)

e = np.array([((0.1*Fc*epc20)/(2*Fc)), 1.25*((0.1*Fc*epc20)/(2*Fc)),
                4.0*((0.1*Fc*epc20)/(2*Fc)), 8.7*((0.1*Fc*epc20)/(2*Fc))])
St = np.array([0.1*Fc, 0.1*0.77*Fc, 0.1*0.45*Fc, 0.1*0.1*Fc])
et = e-((St*epc20)/(2*Fc))  # Applying natural logarithm element-wise

for s, e_val in zip(St, et):
    mct20 +=((round(s, 5), round(e_val, 5)),)

mdb.models[str(name)+'-F-'+str(t6)].Material(name='Concrete-Material-20')
mdb.models[str(name)+'-F-'+str(t6)].materials['Concrete-Material-20'].Density(table=((2.4e-09, ), ))
mdb.models[str(name)+'-F-'+str(t6)].materials['Concrete-Material-20'].Elastic(table=mce20, temperatureDependency=OFF)
mdb.models[str(name)+'-F-'+str(t6)].materials['Concrete-Material-20'].ConcreteDamagedPlasticity(table=((31.0, 0.1, 1.16, 0.667, 0.0001), ))
mdb.models[str(name)+'-F-'+str(t6)].materials['Concrete-Material-20'].concreteDamagedPlasticity.ConcreteCompressionHardening(table=mc20, temperatureDependency=OFF)
mdb.models[str(name)+'-F-'+str(t6)].materials['Concrete-Material-20'].concreteDamagedPlasticity.ConcreteTensionStiffening(table=mct20, temperatureDependency=OFF)

mdb.models[str(name)+'-F-'+str(t6)].sections['Beam'].setValues(material='Concrete-Material-20', thickness=None)
mdb.models[str(name)+'-F-'+str(t6)].sections['Main-Bar'].setValues(area=A3, material='Steel-Material-20')
mdb.models[str(name)+'-F-'+str(t6)].sections['Tie-Bar'].setValues(area=A1, material='Steel-Material-20')
mdb.models[str(name)+'-F-'+str(t6)].sections['Top-Bar'].setValues(area=A2, material='Steel-Material-20')

mdb.models[str(name)+'-F-'+str(t6)].parts['Beam'].setElementType(elemTypes=(
    ElemType(elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT, 
    distortionControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
    ElemType(elemCode=C3D4, elemLibrary=STANDARD)),regions=(
    mdb.models[str(name)+'-F-'+str(t6)].parts['Beam'].cells.findAt(
    ((0.5*b-5, h, S1-5), ), ((0.5*b+5, h, S1-5), ), ((0.5*b-5, 0, S1-5), ), ((0.5*b+5, 0, S1-5), ), 
    ((0.5*b-5, h, S1+5), ), ((0.5*b+5, h, S1+5), ), ((0.5*b-5, 0, S1+5), ), ((0.5*b+5, 0, S1+5), ),
    ((0.5*b-5, h, S2-5), ), ((0.5*b+5, h, S2-5), ), ((0.5*b-5, 0, S2-5), ), ((0.5*b+5, 0, S2-5), ), 
    ((0.5*b-5, h, S2+5), ), ((0.5*b+5, h, S2+5), ), ((0.5*b-5, 0, S2+5), ), ((0.5*b+5, 0, S2+5), ),
    ((0.5*b-5, h, L1-5), ), ((0.5*b+5, h, L1-5), ), ((0.5*b-5, 0, L1-5), ), ((0.5*b+5, 0, L1-5), ), 
    ((0.5*b-5, h, L1+5), ), ((0.5*b+5, h, L1+5), ), ((0.5*b-5, 0, L1+5), ), ((0.5*b+5, 0, L1+5), ),
    ((0.5*b-5, h, L2-5), ), ((0.5*b+5, h, L2-5), ), ((0.5*b-5, 0, L2-5), ), ((0.5*b+5, 0, L2-5), ), 
    ((0.5*b-5, h, L2+5), ), ((0.5*b+5, h, L2+5), ), ((0.5*b-5, 0, L2+5), ), ((0.5*b+5, 0, L2+5), ),
    ), ))

mdb.models[str(name)+'-F-'+str(t6)].parts['Main-Bar'].setElementType(elemTypes=(
    ElemType(elemCode=T3D2, elemLibrary=STANDARD), ), regions=(
    mdb.models[str(name)+'-F-'+str(t6)].parts['Main-Bar'].edges.findAt(((0.5*L, 0.0, 0.0), )), ))

mdb.models[str(name)+'-F-'+str(t6)].parts['Tie-Bar'].setElementType(elemTypes=(
    ElemType(elemCode=T3D2, elemLibrary=STANDARD), ), regions=(
    mdb.models[str(name)+'-F-'+str(t6)].parts['Tie-Bar'].edges.findAt(((0.5*(b-2*Cs-r1), 0.0, 0.0), ), ((0.5*(b-2*Cs-r1), h-Ct-Cb-r1, 0.0), ), 
    ((0.0, 0.5*(h-Ct-Cb-r1), 0.0), ), ((b-2*Cs-r1, 0.5*(h-Ct-Cb-r1), 0.0), ), ), ))

mdb.models[str(name)+'-F-'+str(t6)].parts['Top-Bar'].setElementType(elemTypes=(
    ElemType(elemCode=T3D2, elemLibrary=STANDARD), ), regions=(
    mdb.models[str(name)+'-F-'+str(t6)].parts['Top-Bar'].edges.findAt(((0.5*L, 0.0, 0.0), )), ))

#Job

a = [t6, t1, t2, t3, t4, t5]

for i in range(1, 7):
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model=str(name)+'-F-'+str(a[i-1]), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name='Fire-'+str(a[i-1])+'-min', nodalOutputPrecision=SINGLE, 
            numCpus=8, numDomains=8, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
            '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)


for i in range(1, 7):
        mdb.jobs['Fire-'+str(a[i-1])+'-min'].submit(consistencyChecking=OFF)
        mdb.jobs['Fire-'+str(a[i-1])+'-min'].waitForCompletion()


#lets_Go

from abaqus import *
from abaqusConstants import *
from caeModules import *
import regionToolset 
import visualization
import odbAccess
from driverUtils import executeOnCaeStartup
import os
import numpy as np

times = [t1, t2, t3, t4, t5]
times1 = [t1, t2, t3, t4, t5, t6]

report_path = 'F:\MSc\Post-Fire-' + str(name) + '\'
try:
    os.makedirs(report_path)
    print('Directory created')
except:
    print('Directory exists')

########################post processing#################
for a in times:
	job_name='Fire-'+str(a)+'-min'
	odbname=job_name+'.odb'
	odb = session.openOdb(name=odbname)
	session.viewports['Viewport: 1'].setValues(displayedObject=session.openOdb(name=odbname))
	session.viewports['Viewport: 1'].setValues(displayedObject=session.odbs[odbname])
	session.viewports['Viewport: 1'].setValues(displayedObject=None)
	
	xy1 = xyPlot.XYDataFromHistory(odb=odb, outputVariableName='Reaction force: RF2 PI: rootAssembly Node 3 in NSET RP3', steps=('PF-Load',))
	xy_result1 = session.XYData(name='Reaction-Force', objectToCopy=xy1)
	
	xy2 = xyPlot.XYDataFromHistory(odb=odb, outputVariableName='Spatial displacement: U2 PI: BEAM-1 Node 58 in NSET U', steps=('PF-Load',))
	xy_result2 = session.XYData(name='Displacement', objectToCopy=xy2)
       
	xy3 = session.xyDataObjects['Displacement']
	xy4 = session.xyDataObjects['Reaction-Force']
	xy5 = combine(-xy3, xy4/500)
	xy5.setValues(sourceDescription='')
	tmpName = xy5.name
	
	xy6 = xyPlot.XYDataFromHistory(odb=odb, outputVariableName='Nodal temperature: NT11 PI: MAIN-BAR-1-LIN-1-2 Node 65 in NSET TC-R', steps=('Fire',))
	xy_result1 = session.XYData(name='T-'+str(a)+'-min', objectToCopy=xy6)
	
	session.xyDataObjects.changeKey(tmpName, 'F-'+str(a)+'-min')
    
###write to file

job_name='Fire-0-min'
odbname=job_name+'.odb'
odb = session.openOdb(name=odbname)
session.viewports['Viewport: 1'].setValues(displayedObject=session.openOdb(name=odbname))
session.viewports['Viewport: 1'].setValues(displayedObject=session.odbs[odbname])
session.viewports['Viewport: 1'].setValues(displayedObject=None)
	
xy1 = xyPlot.XYDataFromHistory(odb=odb, outputVariableName='Reaction force: RF2 PI: rootAssembly Node 3 in NSET RP3', steps=('PF-Load',))
xy_result1 = session.XYData(name='Reaction-Force', objectToCopy=xy1)
	
xy2 = xyPlot.XYDataFromHistory(odb=odb, outputVariableName='Spatial displacement: U2 PI: BEAM-1 Node 58 in NSET U', steps=('PF-Load',))
xy_result2 = session.XYData(name='Displacement', objectToCopy=xy2)
       
xy3 = session.xyDataObjects['Displacement']
xy4 = session.xyDataObjects['Reaction-Force']
xy5 = combine(-xy3, xy4/500)
xy5.setValues(sourceDescription='')
tmpName = xy5.name
session.xyDataObjects.changeKey(tmpName, 'F-0-min')

#End

# Modify the report_path variable to the desired folder path
report_path1 = 'C:\Users\DIU\OneDrive - Daffodil International University\MHB-MSc\Post-Fire-' + str(name) + '\Post-Fire-'+ str(name) + '\'

# Ensure the directory exists or create it if it doesn't
try:
    os.makedirs(report_path1)
    print('Directory created')
except:
    print('Directory exists')

# Modify the path when saving CSV files
for b in times1:
    # Retrieve the data from xy5
    data = session.xyDataObjects['F-'+str(b)+'-min'].data
    
    # Save the data to a CSV file
    np.savetxt(report_path1+'-F-'+str(b)+'-min'+'.csv', data, delimiter = ',', fmt = '%.2f', header = 'Displacement, Force')

# Modify the path when saving CSV files
for a in times:
    # Retrieve the data from xy5
    data = session.xyDataObjects['T-'+str(a)+'-min'].data
    
    # Save the data to a CSV file
    np.savetxt(report_path1+'-T-'+str(a)+'-min'+'.csv', data, delimiter = ',', fmt = '%.2f', header = 'Time, Temperature')

#End