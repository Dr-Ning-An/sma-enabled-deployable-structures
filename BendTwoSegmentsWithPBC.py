# -*- coding: mbcs -*-
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
from connectorBehavior import *
#From Python:
import numpy
import math
import datetime


execfile('Parameters_for_BendSingleSegment.py')

L = 40.0;
W = 10.0;
H = t+1.0;
d = 0.1482;

Mdb()

mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=10.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-W/2.0, -H/2.0), 
    point2=(W/2.0, H/2.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseSolidExtrude(depth=L, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
#
mdb.models['Model-1'].ConstrainedSketch(gridSpacing=2.06, name='__profile__', 
    sheetSize=82.48, transform=
    mdb.models['Model-1'].parts['Part-1'].MakeSketchTransform(
    sketchPlane=mdb.models['Model-1'].parts['Part-1'].faces.findAt((-1.666667, 
    -0.166667, L), ), sketchPlaneSide=SIDE1, 
    sketchUpEdge=mdb.models['Model-1'].parts['Part-1'].edges.findAt((W/2.0, 0.0, 
    L), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, L)))
mdb.models['Model-1'].parts['Part-1'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, H/2.0-0.5), point1=(d/2.0, H/2.0-0.5))
mdb.models['Model-1'].parts['Part-1'].PartitionFaceBySketch(faces=
    mdb.models['Model-1'].parts['Part-1'].faces.findAt(((-1.666667, -0.166667, 
    L), )), sketch=mdb.models['Model-1'].sketches['__profile__'], 
    sketchUpEdge=mdb.models['Model-1'].parts['Part-1'].edges.findAt((W/2.0, 0.0, 
    L), ))
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].parts['Part-1'].PartitionCellBySweepEdge(cells=
    mdb.models['Model-1'].parts['Part-1'].cells, 
    edges=(mdb.models['Model-1'].parts['Part-1'].edges.findAt((0.0, 
    H/2.0-0.5+d/2.0, L), ), ), sweepPath=
    mdb.models['Model-1'].parts['Part-1'].edges.findAt((W/2.0, H/2.0, L/2.0), ))

#
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=XYPLANE)
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=XZPLANE)
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=YZPLANE)
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE2, offset=-(H/2.0-0.5)
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[5])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE2, offset=H/2.0-0.5-0.1
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[5])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=0.1
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[8])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=0.2
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[8])
#
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=2.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=5.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=10.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=13.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=18.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=21.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=26.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=29.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=34.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=37.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[4])

#
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=-0.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[6])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=-2.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[6])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=0.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[6])
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByOffset(flip=SIDE1, offset=2.5
    , plane=mdb.models['Model-1'].parts['Part-1'].datums[6])

#
for i in range(5, 25):
    mdb.models['Model-1'].parts['Part-1'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-1'].cells, 
    datumPlane=mdb.models['Model-1'].parts['Part-1'].datums[i])

#
#####
mdb.models['Model-1'].parts['Part-1'].Set(cells=
    mdb.models['Model-1'].parts['Part-1'].cells.getSequenceFromMask(
    mask=('[#0:12 #ffffb800 #ffffff ]', ), ), 
    name='Part-SMA')
mdb.models['Model-1'].parts['Part-1'].Set(cells=
    mdb.models['Model-1'].parts['Part-1'].cells.getSequenceFromMask(mask=(
    '[#2a8055 #1800c #1ffc0fe #30018 #a8000280 #50600032 #c007c8c0', 
    ' #8c00061b #29113 #18000c22 #550022 #c00154 #18 ]', ), ), name='Part-ABS')
mdb.models['Model-1'].parts['Part-1'].Set(cells=
    mdb.models['Model-1'].parts['Part-1'].cells.getSequenceFromMask(mask=(
    '[#ffd57faa #fffe7ff3 #fe003f01 #fffcffe7 #57fffd7f #af9fffcd #3ff8373f', 
    ' #73fff9e4 #fffd6eec #e7fff3dd #ffaaffdd #ff3ffeab #47e7 ]', ), ), name=
    'Part-PDMS')
#
mdb.models['Model-1'].parts['Part-1'].Set(
    edges=mdb.models['Model-1'].parts['Part-1'].edges.getSequenceFromMask(mask=(
    '[#0:46 #408020 #10000000 #800004 #1 #41000 #40000800 ]', ), ), name='Part-SMA-Line')
#: The set 'Part-SMA-Line' has been created (11 edges).
mdb.models['Model-1'].parts['Part-1'].Set(
    edges=mdb.models['Model-1'].parts['Part-1'].edges.getSequenceFromMask(mask=(
    '[#0:17 #4000 #20000004 #20000 #40000200 #400000 #400', ' #1000002 #1 ]', 
    ), ), name='Part-Middle-Line')
#: The set 'Part-Middle-Line' has been created (11 edges).
mdb.models['Model-1'].parts['Part-1'].Set(
    edges=mdb.models['Model-1'].parts['Part-1'].edges.getSequenceFromMask(mask=(
    '[#0:41 #2000000 #1000800 #8000 #4004 #0:4 #800', ' #3000000 #20 #40 ]', ), 
    ), name='Part-ABS-Neutral')
mdb.models['Model-1'].parts['Part-1'].Set(
    faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(mask=(
    '[#0:2 #441418 #a0200041 #20 #0:3 #54000000 #1001010', 
    ' #0:2 #80200800 #2220 #0 #802b0000 #10 #0:2', 
    ' #c0802010 #8 #0 #a8028680 #a000 #22 #0:19', 
    ' #1c8e8000 #6073702e #c3517 ]', ), ), name='Set-OutBound')
mdb.models['Model-1'].parts['Part-1'].Set(
    faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(mask=(
    '[#0:27 #80840000 #60a01 #10802050 #418 #60510212 #8a00000', 
    ' #4181008 #80840000 #60a01 #10802050 #418 #60510212 #8a00000', 
    ' #4181008 #40840000 #50a01 #40804c0 #20100428 #f080010 ]', ), ), name='Fix_Z')
######
#
mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.rotate(instanceList=('Part-1-1', ), axisPoint=(0.0, 0.0, 0.0), 
    axisDirection=(0.0, 1.0, 0.0), angle=90.0)
mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
    0.0), direction2=(0.0, 1.0, 0.0), instanceList=('Part-1-1', ), number1=2, 
    number2=1, spacing1=L, spacing2=3.0)
mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(1.0, 0.0, 
    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-1-1-lin-2-1', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-1-lin-2-1', 
    ), vector=(0.0, t, 0.0))
#
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1-lin-2-1']), 
    keepIntersections=ON, name='Part-2', originalInstances=SUPPRESS)
#
mdb.models['Model-1'].rootAssembly.makeIndependent(instances=(mdb.models['Model-1'].rootAssembly.instances['Part-2-1'], ))

#
#
mdb.models['Model-1'].rootAssembly.Set(name='Fix_X', vertices=
    mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].vertices.findAt(((
    0.0, 0.0, 0.0), )))
##
#
##
if H==1.5:
    mdb.models['Model-1'].rootAssembly.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
        plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=0.5)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
        plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=0.65)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
        plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=-0.15)
    for i in range(8,12):
        mdb.models['Model-1'].rootAssembly.PartitionCellByDatumPlane(
        datumPlane=mdb.models['Model-1'].rootAssembly.datums[i], 
        cells=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells)
#
if H==1.6:
    mdb.models['Model-1'].rootAssembly.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
        plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=t)
    for i in range(8,10):
        mdb.models['Model-1'].rootAssembly.PartitionCellByDatumPlane(
        datumPlane=mdb.models['Model-1'].rootAssembly.datums[i], 
        cells=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells)
#
if H>1.6 and H<2.0:
    mdb.models['Model-1'].rootAssembly.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
        plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=t)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
        plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=H/2.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
        plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=-H/2.0+t)
    for i in range(8,12):
        mdb.models['Model-1'].rootAssembly.PartitionCellByDatumPlane(
        datumPlane=mdb.models['Model-1'].rootAssembly.datums[i], 
        cells=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells)
if H>2.0:
#
    mdb.models['Model-1'].rootAssembly.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
    	plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=(H-1.0-1.0)/2.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
    	plane=mdb.models['Model-1'].rootAssembly.datums[9], flip=SIDE1, offset=1.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
    	plane=mdb.models['Model-1'].rootAssembly.datums[10], flip=SIDE1, offset=H-1.0-1.0)
    mdb.models['Model-1'].rootAssembly.DatumPlaneByOffset(
		plane=mdb.models['Model-1'].rootAssembly.datums[8], flip=SIDE1, offset=-(H-1.0-1.0)/2.0)
    mdb.models['Model-1'].rootAssembly.PartitionCellByDatumPlane(
        datumPlane=mdb.models['Model-1'].rootAssembly.datums[9], 
        cells=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells)
    mdb.models['Model-1'].rootAssembly.PartitionCellByDatumPlane(
        datumPlane=mdb.models['Model-1'].rootAssembly.datums[10], 
        cells=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells)
    mdb.models['Model-1'].rootAssembly.PartitionCellByDatumPlane(
        datumPlane=mdb.models['Model-1'].rootAssembly.datums[11], 
        cells=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells)
    mdb.models['Model-1'].rootAssembly.PartitionCellByDatumPlane(
        datumPlane=mdb.models['Model-1'].rootAssembly.datums[12], 
        cells=mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells)
#
#
mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
    minSizeFactor=0.1, regions=(
    mdb.models['Model-1'].rootAssembly.instances['Part-2-1'], ), size=0.5)
#
#
mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
    mdb.models['Model-1'].rootAssembly.instances['Part-2-1'], ))
mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    kinematicSplit=AVERAGE_STRAIN, hourglassControl=ENHANCED, 
    distortionControl=DEFAULT),), regions=(
    mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].cells, ))

#
elementsAll_SMA=mdb.models['Model-1'].rootAssembly.sets['Part-2-1.Part-SMA'].elements

# Material
mdb.models['Model-1'].Material(name='Material-SMA')
mdb.models['Model-1'].materials['Material-SMA'].Depvar(n=100)
mdb.models['Model-1'].materials['Material-SMA'].UserMaterial(mechanicalConstants=
    (2.0, 3.0, 1e-20, 1.0, len(elementsAll_SMA), 75000.0, 28000.0, 0.33, 2.2e-05, 2.2e-05, 
    350.0, 315.0, 354.0, 379.0, 0.065, -0.35, -0.35, 0.07, -0.035, -0.035, 0.0, 0.0, 
    0.0, 1.0))
mdb.models['Model-1'].materials['Material-SMA'].Density(table=((6.45e-9, ), 
    ))
mdb.models['Model-1'].HomogeneousSolidSection(material='Material-SMA', name=
    'Section-SMA', thickness=None)
mdb.models['Model-1'].parts['Part-2'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-2'].sets['Part-SMA'], sectionName=
    'Section-SMA', thicknessAssignment=FROM_SECTION)
#
mdb.models['Model-1'].Material(name='Material-PDMS')
mdb.models['Model-1'].materials['Material-PDMS'].Elastic(table=((1.8, 0.495), ))
mdb.models['Model-1'].materials['Material-PDMS'].Density(table=((0.965e-9, ), 
    ))
mdb.models['Model-1'].HomogeneousSolidSection(material='Material-PDMS', name=
    'Section-PDMS', thickness=None)
#
mdb.models['Model-1'].parts['Part-2'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-2'].sets['Part-PDMS'], sectionName=
    'Section-PDMS', thicknessAssignment=FROM_SECTION)
#
mdb.models['Model-1'].Material(name='Material-ABS')
mdb.models['Model-1'].materials['Material-ABS'].Elastic(table=((2.6e3, 0.35), ))
mdb.models['Model-1'].materials['Material-ABS'].Density(table=((1.2e-9, ), 
    ))
mdb.models['Model-1'].HomogeneousSolidSection(material='Material-ABS', name=
    'Material-ABS', thickness=None)
mdb.models['Model-1'].parts['Part-2'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-2'].sets['Part-ABS'], sectionName=
    'Material-ABS', thicknessAssignment=FROM_SECTION)

#
mdb.models['Model-1'].ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, 
    application=QUASI_STATIC, initialConditions=OFF, initialInc=0.005, maxInc=
    0.005, maxNumInc=1000, name='Step-1', nlgeom=ON, nohaf=OFF, previous=
    'Initial')
mdb.models['Model-1'].FieldOutputRequest(name='F-Output-1', 
    createStepName='Step-1', variables=PRESELECT)
mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-1', 
    createStepName='Step-1', variables=PRESELECT)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'V', 'A', 'RF', 'CF', 'CSTRESS', 
    'CDISP', 'COORD', 'SDV'))
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(
    frequency=1)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
    frequency=1)
#
NameModel='Model-1'; NameRef1='RefPoint-0'; NameSet='Part-2-1.Set-OutBound'
mdb.models[NameModel].Part(dimensionality=THREE_D, name=NameRef1, type=
    DEFORMABLE_BODY)
mdb.models[NameModel].parts[NameRef1].ReferencePoint(point=(0.0, 0.0, 0.0))
mdb.models[NameModel].rootAssembly.Instance(dependent=ON, name=NameRef1, 
    part=mdb.models[NameModel].parts[NameRef1])
#Create set of reference points
mdb.models[NameModel].rootAssembly.Set(name=NameRef1, referencePoints=(
    mdb.models[NameModel].rootAssembly.instances[NameRef1].referencePoints[1],))
nodesAll=mdb.models[NameModel].rootAssembly.sets[NameSet].nodes
nodesAllCoor=[]
for nod in mdb.models[NameModel].rootAssembly.sets[NameSet].nodes:
    nodesAllCoor.append(nod.coordinates)
repConst=0
#Find periodically located nodes and apply equation constraints
ranNodes=range(0,len(nodesAll)) #Index array of nodes not used in equations constraint
for repnod1 in range(0,len(nodesAll)):
    stop=False          #Stop will become true when equation constraint is made between nodes
    Coor1=nodesAllCoor[repnod1]     #Coordinates of Node 1
    for repnod2 in ranNodes:    #Loop over all available nodes
        Coor2=nodesAllCoor[repnod2] #Coordinates of Node 2
        dx=Coor2[0]-Coor1[0]; dy=Coor2[1]-Coor1[1]; dz=Coor2[2]-Coor1[2]    #X and Y Distance between nodes
        if int(round(1000.0*dy))==0 and int(round(1000.0*dz))==0 and int(round(1000.0*(2.0*L-dx)))==0:
            mdb.models[NameModel].rootAssembly.Set(name='Node-1-'+str(repConst), nodes=
                mdb.models[NameModel].rootAssembly.sets[NameSet].nodes[repnod1:repnod1+1])
            mdb.models[NameModel].rootAssembly.Set(name='Node-2-'+str(repConst), nodes=
                mdb.models[NameModel].rootAssembly.sets[NameSet].nodes[repnod2:repnod2+1])
            for Dim1 in [1,2]:
                mdb.models[NameModel].Equation(name='PerConst'+str(Dim1)+'-'+str(repConst),
                    terms=((1.0,'Node-1-'+str(repConst), Dim1),(-1.0, 'Node-2-'+str(repConst), Dim1) ,
                    (1.0, 'RefPoint-0', Dim1)))
            repConst=repConst+1 #Increase integer for naming equation constraint
            ranNodes.remove(repnod1)#Remove used node from available list
            stop=True       #Don't look further, go to following node.
#
#
mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
    region=mdb.models['Model-1'].rootAssembly.sets['RefPoint-0'], u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
    region=mdb.models['Model-1'].rootAssembly.sets['Fix_X'], u1=SET, 
    u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-3', 
    region=mdb.models['Model-1'].rootAssembly.sets['Part-2-1.Fix_Z'], u1=UNSET, 
    u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
##
mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(353.0, ), name='Predefined Field-1', region=
    mdb.models['Model-1'].rootAssembly.sets['Part-2-1.Part-SMA'])
mdb.models['Model-1'].predefinedFields['Predefined Field-1'].setValuesInStep(
    magnitudes=(420.0, ), stepName='Step-1')

# Job
# os.chdir('L%dt%d'%(round(L), round(10*t)))
curr = os.getcwd()
JobName = 'BendTwoSegmentsWithPBC'
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=JobName, nodalOutputPrecision=SINGLE, 
    numCpus=2, numDomains=2, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine=curr + '\\sma_um.for'
    , waitHours=0, waitMinutes=0)
mdb.saveAs(pathName=JobName)
mdb.jobs[JobName].submit(consistencyChecking=OFF)
mdb.jobs[JobName].waitForCompletion()

# PostProcessing

from odbAccess import*
from abaqusConstants import * 
import string   
import numpy as np
import os

odb = openOdb(path = JobName+'.odb')

os.makedirs('ResultFiles')
os.chdir('ResultFiles')

for fm in range(0, len(odb.steps['Step-1'].frames)):
    outfile = open('SMA-X-Y-Z-Temp-' + str(fm) + '.csv','w')
    outfile.write('X, Y, Z, Temperature\n')
    timeFrame = odb.steps['Step-1'].frames[fm]
    readNode = odb.rootAssembly.instances['PART-2-1'].nodeSets['PART-SMA-LINE']
    Coordinate = timeFrame.fieldOutputs['COORD']
    readNodeCoordinate = Coordinate.getSubset(region=readNode)
    readNodeCoordinateValues = readNodeCoordinate.values
    count=len(readNodeCoordinateValues)
    X_Coordinate = np.zeros(count)
    Y_Coordinate = np.zeros(count)
    Z_Coordinate = np.zeros(count)
    for i in range(0, count):
        X_Coordinate[i]=readNodeCoordinateValues[i].data[0]
        Y_Coordinate[i]=readNodeCoordinateValues[i].data[1]
        Z_Coordinate[i]=readNodeCoordinateValues[i].data[2]
    Sorted_X_Coordinate = np.sort(X_Coordinate)
    Inps = X_Coordinate.argsort()
    Sorted_Y_Coordinate = Y_Coordinate[Inps]
    Sorted_Z_Coordinate = Z_Coordinate[Inps]
    for i in range(0, count):
        outfile.write(str(Sorted_X_Coordinate[i]) + ',' +
        str(Sorted_Y_Coordinate[i]) + ',' + 
        str(Sorted_Z_Coordinate[i]) + ',' +'\n')
    outfile.close()

for fm in range(0, len(odb.steps['Step-1'].frames)):
    outfile = open('ABS-X-Y-Z-Temp-' + str(fm) + '.csv','w')
    outfile.write('X, Y, Z, Temperature\n')
    timeFrame = odb.steps['Step-1'].frames[fm]
    readNode = odb.rootAssembly.instances['PART-2-1'].nodeSets['PART-ABS-NEUTRAL']
    Coordinate = timeFrame.fieldOutputs['COORD']
    readNodeCoordinate = Coordinate.getSubset(region=readNode)
    readNodeCoordinateValues = readNodeCoordinate.values
    count=len(readNodeCoordinateValues)
    X_Coordinate = np.zeros(count)
    Y_Coordinate = np.zeros(count)
    Z_Coordinate = np.zeros(count)
    for i in range(0, count):
        X_Coordinate[i]=readNodeCoordinateValues[i].data[0]
        Y_Coordinate[i]=readNodeCoordinateValues[i].data[1]
        Z_Coordinate[i]=readNodeCoordinateValues[i].data[2]
    Sorted_X_Coordinate = np.sort(X_Coordinate)
    Inps = X_Coordinate.argsort()
    Sorted_Y_Coordinate = Y_Coordinate[Inps]
    Sorted_Z_Coordinate = Z_Coordinate[Inps]
    for i in range(0, count):
        outfile.write(str(Sorted_X_Coordinate[i]) + ',' +
        str(Sorted_Y_Coordinate[i]) + ',' + 
        str(Sorted_Z_Coordinate[i]) + ',' +'\n')
    outfile.close()

for fm in range(0, len(odb.steps['Step-1'].frames)):
    outfile = open('Middle-X-Y-Z-Temp-' + str(fm) + '.csv','w')
    outfile.write('X, Y, Z, Temperature\n')
    timeFrame = odb.steps['Step-1'].frames[fm]
    readNode = odb.rootAssembly.instances['PART-2-1'].nodeSets['PART-MIDDLE-LINE']
    Coordinate = timeFrame.fieldOutputs['COORD']
    readNodeCoordinate = Coordinate.getSubset(region=readNode)
    readNodeCoordinateValues = readNodeCoordinate.values
    count=len(readNodeCoordinateValues)
    X_Coordinate = np.zeros(count)
    Y_Coordinate = np.zeros(count)
    Z_Coordinate = np.zeros(count)
    for i in range(0, count):
        X_Coordinate[i]=readNodeCoordinateValues[i].data[0]
        Y_Coordinate[i]=readNodeCoordinateValues[i].data[1]
        Z_Coordinate[i]=readNodeCoordinateValues[i].data[2]
    Sorted_X_Coordinate = np.sort(X_Coordinate)
    Inps = X_Coordinate.argsort()
    Sorted_Y_Coordinate = Y_Coordinate[Inps]
    Sorted_Z_Coordinate = Z_Coordinate[Inps]
    for i in range(0, count):
        outfile.write(str(Sorted_X_Coordinate[i]) + ',' +
        str(Sorted_Y_Coordinate[i]) + ',' + 
        str(Sorted_Z_Coordinate[i]) + ',' +'\n')
    outfile.close()

odb.close()