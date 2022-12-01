import bpy
import numpy as np
import os 
os.system('cls')#clear console


coll=bpy.data.collections['Collection'] #coll is a variable that encompasses all the objects in a scene

for ob in coll.objects:
    bpy.data.objects.remove(ob)#deletes all objects
    

datos_bichos=np.loadtxt('D:\\un\\metodos\\proyecto\\bichines\\Bichin.txt')
datos_comida=np.loadtxt('D:\\un\\metodos\\proyecto\\bichines\\Food.txt')

#datos has the information array that we need. The first column is time
#Where every row is a certain time, every other ther row is part of an xy pair
#hence you should always have an odd number of columns in datos

sp=datos_bichos.shape
num_obj=(sp[1]-1)/3 # to account for 1 time column 

t=datos_bichos[:,0]

print(f'Tenemos {len(t)} pasos de tiempo con {num_obj} bichines')

def ipos(i,axis): #where i is the object number you want, N => (1,num_obj)
    if(axis=='x'):
        index=3*i-2
    if(axis=='y'):
        index=3*i-1
    if(axis=='l'):
        index=3*i
    return index
    
#-------creation of the objects-----------#

for i in range(1,int(num_obj+1)):
    s=datos_bichos[0,ipos(i,'l')]
    bpy.ops.mesh.primitive_uv_sphere_add(radius=1, location=( datos_bichos[0,ipos(i,'x')]/5, datos_bichos[0,ipos(i,'y')]/5, 0), scale=(s, s, s))
    
    
    mat = bpy.data.materials.get("Material.002")
    if mat is None:
        # create material
        mat = bpy.data.materials.new(name="Material.002")
        bpy.data.materials["Material.002"].node_tree.nodes["Principled BSDF"].inputs[0].default_value = (0.8, 0.122186, 0.0280791, 1)
    ob = bpy.context.active_object
    if ob.data.materials:
    # assign to the first material slot
        ob.data.materials[0] = mat
    else:
        # no slots
        ob.data.materials.append(mat)

    bpy.ops.object.shade_smooth()
    bpy.ops.object.hide_viewport

    #this initializaces a series of spheres, one for each column triad.
    

coll=bpy.data.collections['Collection']
#print(datos_bichos)

#-------------animation---------------#

for time in range(len(t)//2):
    i=0
    for ob in coll.objects:
        i+=1
        #print(datos_bichos[time,ipos(i,'l')])
        if datos_bichos[time,ipos(i,'l')]==1: #if it is alive
            s=datos_bichos[time,ipos(i,'l')]
            bpy.context.scene.objects[ob.name].location=(datos_bichos[time,ipos(i,'x')]/5,datos_bichos[time,ipos(i,'y')]/5,0)
            bpy.context.scene.objects[ob.name].scale=(s,s,s)
            bpy.context.scene.objects[ob.name].keyframe_insert(data_path = 'scale',  frame=float(time*4))
            bpy.context.scene.objects[ob.name].keyframe_insert(data_path='location', frame=float(time*4))
            
            if datos_bichos[time-1,ipos(i,'l')]==0: #if it was just born
                print(f"un nacimiento, t_a={datos_bichos[time-1,ipos(i,'l')]} y t_b={datos_bichos[time,ipos(i,'l')]}  en tl tiempo {time} ")
                
                bpy.context.scene.objects[ob.name].location=(datos_bichos[time,ipos(i,'x')]/5,datos_bichos[time,ipos(i,'y')]/5,0)
                bpy.context.scene.objects[ob.name].scale=(0,0,0)
                bpy.context.scene.objects[ob.name].keyframe_insert(data_path = 'scale',  frame=float(time*4-1))
                bpy.context.scene.objects[ob.name].keyframe_insert(data_path='location', frame=float(time*4-1))
        
