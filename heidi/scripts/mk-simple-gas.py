#!/usr/bin/python

num_bodies = 1000
num_types = 2

def print_particle(type, Z, N, E, x,y,z, vx,vy,vz, ax,ay,az):
    print "%i %i %i %i %e %e %e  %e %e %e  %e %e %e" % (type, Z, N, E, x,y,z, vx,vy,vz, ax,ay,az)

print "NumBodies:", num_bodies
print "NumTypes:",  num_types
print

bounds = (-1

print_particle(1, 79, 79, 0, 0,0,0, 0,0,0, 0,0,0)

start_y = (num_bodies-2) / 2 * 1e24
for i in range(0, num_bodies-1):
    print_particle(2, 2, 2, 0, -2e36 - i*1e36,start_y - i * 1e24,0, 1.03e-01,0,0, 0,0,0)

print

