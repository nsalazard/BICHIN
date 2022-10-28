import bpy
import numpy as np
import os 
os.system('cls')#clear console


coll=bpy.data.collections['Collection'] #coll is a variable that encompasses all the objects in a scene

for ob in coll.objects:
    bpy.data.objects.remove(ob)#deletes all objects
    

datos=np.loadtxt('D:\\un\\metodos\\proyecto\\pruebas_blend\\datos.dat')

#datos has the information array that we need. The first column is time
#Where every row is a certain time, every other ther row is part of an xy pair
#hence you should always have an odd number of columns in datos

sp=datos.shape
num_obj=(sp[1]-1)/2 # to account for 1 time column 

t=datos[:,0]



def ipos(i,axis): #where i is the object number you want, N => (1,num_obj)
    if(axis=='x'):
        index=2*i-1
    if(axis=='y'):
        index=2*i
    return index
    
#-------creation of the objects-----------#

for i in range(1,int(num_obj+1)):
    bpy.ops.mesh.primitive_uv_sphere_add(radius=1, location=( datos[0,ipos(i,'x')]/10, datos[0,ipos(i,'y')]/10, 0), scale=(1, 1, 1))
    #this initializaces a series of spheres, one for each column pair.
    

coll=bpy.data.collections['Collection']

#-------------animation---------------#

for time in range(500):
    bpy.context.scene.frame_set(time)
    i=0
    for ob in coll.objects:
        i+=1
        bpy.context.scene.objects[ob.name].location=(datos[time,ipos(i,'x')]/10,datos[time,ipos(i,'y')]/10,0)
        bpy.context.scene.objects[ob.name].keyframe_insert(data_path='location', frame=float(time))
        
        
