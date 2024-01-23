#IMPORTS OF OUTSIDE LIBRARIES----------

import numpy as np
import matplotlib.pyplot as plt
import random 
import pygame
import pymunk
import pymunk.pygame_util
from pymunk.vec2d import Vec2d
import math
#PYGAME COLORS-------------------------
white =     ( 255, 255, 255 )
blue =      ( 0, 0, 255 )
aliceblue=  (240,248,255)
cobalt=     (61,89,171)
green =     ( 0, 255, 0 )
red =       ( 255, 0, 0 )
black=      ( 0, 0, 0 )
gray17=     ( 43, 43, 43)
#-------------
#----------------------
# SPACE:
space=pymunk.Space()
space.damping=0.75
# PARTICLES:
#cherry 
cherry = pymunk.Body()
cherry.position=400,700
shape= pymunk.Circle(cherry, 13)
shape.mass = 3
shape.elasticity = 0.0
cherry.angular_velocity=40
shape.friction = 0.9 
space.add(cherry, shape)
#-------
#particle1
body=pymunk.Body(body_type=pymunk.Body.DYNAMIC)
body.position=400,400
bshape=pymunk.Circle(body,15)
bshape.mass=100
bshape.elasticity=0
bshape.friction=0.5
space.add(body,bshape)
#--------------
#particle2
part=pymunk.Body()
part.position=415,415
pshape=pymunk.Circle(part,13)
pshape.mass=3
pshape.elasticity=0.7
pshape.friction=0.5
space.add(part,pshape)
#-----------------
#particle3 (proton)
pro=pymunk.Body()
pro.position=(385,385)
proshape=pymunk.Circle(pro,13)
proshape.mass=5
proshape.elasticity=0.7
proshape.friction=0.5
space.add(pro,proshape)
#-----------------


# DEFINED FUNCTIONS:

#Window Boundaries


def boundary_set(space, width, height): #same output of code, different format for more precise addition and change of code
    rects = [
        [(width/2, height ), (width, 10)],
        [(width/2, 0), (width, 10)],
        [(0, height/2), (10, height)],
        [(width , height/2), (10, height)]
    ]
    
    for pos, size in rects: #we create a "body" for the border that is static so we can also add things to the border that may affect the particles state (such as the elasticity and friction of it)
        body = pymunk.Body(body_type=pymunk.Body.STATIC)
        body.position = pos
        shape = pymunk.Poly.create_box(body,size)
        shape.elasticity = 0.5  #less elasticity on border will result in a drop off in "bounciness"
        shape.friction = 0.5
        shape.color=(230,230,250,100)
        space.add(body, shape)  
#----------------------

#Magnetic Forces:

def magnetic_attractive_force(body1,body2,p1,p2,charge1,charge2)-> tuple:
    # Distance (r):
    r=(math.sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2))
    
    # Magnitude of Force:
    k=9*10**9 #Nm^2/c^2
    f=k*(charge1*charge2)/r**2
   
    #----------
    
    # Direction of Force:

    #if p1 is at the same y-direction as p2
    if p2[1]==p1[1]: 
        theta=0
        if p2[0]>p1[0]:
            F1,F2 = (f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        
        elif p2[0]<p1[0]:
            F1,F2 = (-f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    #if p1 is at the same x-direction as p2
    elif p2[0]==p1[0]: 
        theta=np.pi/2
        if p2[1]>p1[1]:
            F1,F2 = (f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        elif p2[1]<p1[1]:
            F1,F2 = (f*np.cos(theta)),(-f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    #--------
    #p1 IS ABOVE p2 
    elif p1[1]>p2[1]:   
        if p2[0]>p1[0]: #if p1 is to the left of p2
            theta=(np.arctan((p2[0]-p1[0])/(p1[1]-p2[1])))
            F1,F2 = (f*np.cos(theta)),(-f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        elif p2[0]<p1[0]: #if p1 is to the right of p2
            theta=(np.arctan((p1[0]-p2[0])/(p1[1]-p2[1])))
            F1,F2 = (-f*np.cos(theta)),(-f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    # p1 IS BELOW p2
    elif p1[1]<p2[1]: 
        if p2[0]>p1[0]: #if p1 is to the left of p2
            theta=(np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])))
            F1,F2 = (f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        elif p2[0]<p1[0]: # if p1 is to the right of p2
            theta=(np.arctan((p2[1]-p1[1])/(p1[0]-p2[0])))
            F1,F2 = (-f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    
    return()
def magnetic_repulsive_force(body1,body2,p1,p2,charge1,charge2)-> tuple:
    # Distance (r):
    r=(math.sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2))
    
    # Magnitude of Force:
    k=9*10**9 #Nm^2/c^2
    f=k*(charge1*charge2)/r**2
   
    #----------
    
    # Direction of Force:
    #if p1 is at the same y-direction as p2
    if p2[1]==p1[1]: 
        theta=0
        if p2[0]>p1[0]:
            F1,F2 = (-f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        elif p2[0]<p1[0]:
            F1,F2 = (f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    #if p1 is at the same x-direction as p2
    elif p2[0]==p1[0]: 
        theta=np.pi/2
        if p2[1]>p1[1]:
            F1,F2 = (f*np.cos(theta)),(-f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        elif p2[1]<p1[1]:
            F1,F2 = (f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    #--------
    #p1 IS ABOVE p2
    elif p1[1]>p2[1]:    
        if p2[0]>p1[0]: #if p1 is to the left of p2
            theta=(np.arctan((p2[0]-p1[0])/(p1[1]-p2[1])))
            F1,F2 = (-f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        elif p2[0]<p1[0]: #if p1 is to the right of p2
            theta=(np.arctan((p1[0]-p2[0])/(p1[1]-p2[1])))
            F1,F2 = (f*np.cos(theta)),(f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    # p1 IS BELOW p2
    elif p1[1]<p2[1]: 
        if p2[0]>p1[0]: #if p1 is to the left of p2
            theta=(np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])))
            F1,F2 = (-f*np.cos(theta)),(-f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
        elif p2[0]<p1[0]: # if p1 is to the right of p2
            theta=(np.arctan((p2[1]-p1[1])/(p1[0]-p2[0])))
            F1,F2 = (f*np.cos(theta)),(-f*np.sin(theta))
            if body2==None:
                return (body1.apply_force_at_world_point((F1,F2),(p2)))
            else:
                return (body1.apply_force_at_world_point((F1,F2),(p2)),body2.apply_force_at_world_point((-F1,-F2),(p1)))
    
    return()  
#--------------------------------------------


pygame.init()
WIDTH, HEIGHT = 800, 800
window = pygame.display.set_mode((WIDTH, HEIGHT))

def run (window, width, height): #just redid our main loop 
    run = True
    clock = pygame.time.Clock()
    fps = 60
    dt = 1/fps

    while run:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False
                break
        window.fill((gray17))

        magnetic_repulsive_force(cherry,part,cherry.position,part.position,0.01,0.01)
        magnetic_attractive_force(cherry,body,cherry.position,body.position,0.02,0.02)
        magnetic_attractive_force(part,body,part.position,body.position,0.015,0.015)
    



       
        pygame.draw.circle(window,aliceblue,(cherry.position),13,0)
        pygame.draw.line(window,blue,(cherry.position[0]-7,cherry.position[1]),(cherry.position[0]+7,cherry.position[1]),3)
        pygame.draw.circle(window,aliceblue,(part.position),13,0)
        pygame.draw.line(window,blue,(part.position[0]-7,part.position[1]),(part.position[0]+7,part.position[1]),3)
        pygame.draw.circle(window,cobalt,(body.position),15,0)
        pygame.draw.line(window,white,(body.position[0]-7,body.position[1]),(body.position[0]+7,body.position[1]),3)
        pygame.draw.line(window,white,(body.position[0],body.position[1]-7),(body.position[0],body.position[1]+7),3)
        
        
        #print("KE_cherry:",cherry.kinetic_energy)
    
    
        space.step(dt)
        clock.tick(fps)
        boundary_set(space, width, height)
        pygame.display.update()
    pygame.quit()
    
if __name__ == "__main__":
    run(window, WIDTH, HEIGHT)
