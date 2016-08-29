#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Name:        tubewithcaps.py
# Purpose:     Generate tubes in blender which follow a path. The tubes have a
#              circular profile, and hemispherical end caps.
#
# Author:      Brian G. Olson (olsonbg@gmail.com)
#
# Created:     27 October 2015
# Copyright:   (c) 2015
# -----------------------------------------------------------------------------

import bpy
import bmesh
import mathutils


class TubeWithCaps(object):

    def __init__(self, name, radius, coordinates):
        self.name         = name
        self.radius       = radius
        self.coordinates  = coordinates
        self.startCapName = "CapStart"
        self.endCapName   = "CapEnd"
        self.pathName     = "Path"
        self.bevelName    = "ChainBevel"

    #
    # Bevel object for curves
    #
    def make_bevelobject(self):
        # Circle for bevel object of curves (tubes).
        bpy.ops.curve.primitive_bezier_circle_add(radius=self.radius,
                                                  location=(0, 0, 0))
        bpy.context.object.name = self.bevelName

        return bpy.data.objects[self.bevelName]

    # Hemisphere on positive z.
    def makeHemisphere(self, objname, segments):

        ob = bpy.ops.mesh.primitive_uv_sphere_add(segments=segments,
                                                  location=(0, 0, 0),
                                                  size=self.radius)
        bpy.ops.object.shade_smooth()
        bpy.context.object.name = objname

        bpy.ops.object.mode_set(mode='EDIT')
        mesh = bmesh.from_edit_mesh(bpy.context.object.data)
        for v in mesh.verts:
            if v.co[2] <  0.0:
                v.select = True
            else:
                v.select = False
                #  mesh.verts.remove(v)

        bpy.ops.mesh.delete(type='VERT')
        bpy.ops.object.mode_set(mode='OBJECT')

        return bpy.data.objects[objname]

    def make_caps(self, num_hemisphere_verts):
        segments = num_hemisphere_verts
        startCap = self.makeHemisphere(self.startCapName, segments)

        # Make a copy for the endCap
        endCap          = startCap.copy()
        endCap.location = startCap.location
        endCap.data     = startCap.data.copy()  # duplicate mesh
        endCap.name     = self.endCapName

        bpy.context.scene.objects.link(endCap)

        return startCap, endCap

    def merge_verts(self, tube):
        bm = bmesh.new()

        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(tube.data)  # Fill bmesh with object mesh

        # get vertices of merged tube and caps which have 3 connections. They
        # should be on the ends of the tube, and the caps
        tube_end_verts = self.verts_withNconnections(bm, 3)

        num_verts = int(len(tube_end_verts)/4)

        v1 = tube_end_verts[:num_verts]
        v2 = tube_end_verts[3*num_verts:4*num_verts]
        v3 = tube_end_verts[num_verts:2*num_verts]
        v4 = tube_end_verts[2*num_verts:3*num_verts]

        # v1, v2, v3, and v4 should be same length
        for i in range(len(v1)):
            # v1 and v2 run in opposite directions, so use j for v2 and i for
            # v1.
            j = len(v1)-i-1+1
            if i == 0:
                j = 0

            # Put the vertices which should be merged at the same location so
            # that remove_doubles will have no problem finding them
            v2[j].co = v1[i].co
            v4[i].co = v3[i].co

        bmesh.ops.remove_doubles(bm, verts=v1+v2)
        bmesh.ops.remove_doubles(bm, verts=v3+v4)
        bpy.ops.object.mode_set(mode='OBJECT')

    # Get all vertices which have N edges.
    def verts_withNconnections(self, bm, N):

        verts = [v for v in bm.verts if len(v.link_edges) == N]

        return verts

    def verts_withNconnections_co(self, ob, N):
        # Use bmesh
        bm = bmesh.new()
        bm.from_mesh(ob.data)  # Fill bmesh with object mesh

        # Get all vertices which have N edges.
        verts = self.verts_withNconnections(bm, N)

        # Get the coordinates of the the vertices
        verts_co = [(ob.matrix_world * v.co) for v in verts]

        return verts_co

    def extrude_tube(self):

        polyline, curvedata = self.make_polyline()

        bevel_ob = self.make_bevelobject()
        curvedata.bevel_object = bevel_ob

        # Convert to a mesh
        bpy.data.objects[self.name].select = True
        bpy.ops.object.convert(target='MESH')

        # Done with bevel object, so delete it.
        bpy.ops.object.select_all(action='DESELECT')
        bevel_ob.select = True
        bpy.ops.object.delete()

        return polyline
        #  return bpy.data.objects[self.name]

    def make_polyline(self):

        curvedata = bpy.data.curves.new(name=self.pathName, type='CURVE')
        curvedata.dimensions = '3D'

        objectdata = bpy.data.objects.new(self.name, curvedata)
        objectdata.location = (0, 0, 0)  # object origin
        bpy.context.scene.objects.link(objectdata)

        polyline = curvedata.splines.new('NURBS')
        polyline.points.add(len(self.coordinates)-1)

        w = 100.0  # Weight of control points
        for i, coords in enumerate(self.coordinates):
            x, y, z = coords
            polyline.points[i].co = (x, y, z, w)

        polyline.order_u = len(polyline.points)-1

        # Make curve touch endpoints
        polyline.use_endpoint_u = True

        return objectdata, curvedata

    def aligncap(self, tube_end_co, cap, cap_normal, location):

        # Get normal to tube end. Pick three points equally spaced to determine
        # normals
        one_third = int(len(tube_end_co)/3)
        end = [tube_end_co[0],
               tube_end_co[one_third],
               tube_end_co[2*one_third]]

        normal = mathutils.geometry.normal(*end)
        # normal_sphere = mathutils.Vector((0, 0, -1))
        to_rotate = cap_normal.rotation_difference(normal)
        cap.rotation_mode = 'QUATERNION'
        # Rotate the cap so that it is at the same angle as the tube ends
        cap.rotation_quaternion = to_rotate

        # Move cap to the tube end
        cap.location = location

        # Apply transforms for later calculations.
        bpy.ops.object.transform_apply(scale=True)

        # Rotate cap to align vertices with tube vertices. The vertices will
        # be very close to aligned, but not perfectly.
        capverts = self.verts_withNconnections_co(cap, 3)
        v1 = capverts[0]    - mathutils.Vector(location)
        v2 = tube_end_co[0] - mathutils.Vector(location)
        anglediff = v1.rotation_difference(v2)
        cap.delta_rotation_quaternion = -1.0*anglediff

    def make(self):
        #  bpy.data.objects["ChainBevel"].select = False
        tube = self.extrude_tube()
        #  tube = bpy.data.objects[self.name]

        # Get vertices of tube which have 3 connections. They should be on the
        # ends of the tube.
        tube_end_verts_co = self.verts_withNconnections_co(tube, 3)
        num_verts = int(len(tube_end_verts_co)/2)

        # Split verts in half. The vStart, and second half vEnd
        tube_start_co = tube_end_verts_co[:num_verts]
        tube_end_co   = tube_end_verts_co[num_verts:]

        # Make hemispheres with same number of segments as the tube
        capStart, capEnd = self.make_caps(num_verts)

        self.aligncap(tube_start_co,
                      capStart,
                      mathutils.Vector((0, 0, -1)),
                      self.coordinates[0])

        self.aligncap(tube_end_co,
                      capEnd,
                      mathutils.Vector((0, 0, 1)),
                      self.coordinates[-1])

        # Join caps and tube
        # The active object will be the name of the joined objects.
        capEnd.select   = True
        tube.select     = True
        capStart.select = True
        bpy.context.scene.objects.active = tube

        bpy.ops.object.join()

        self.merge_verts(tube)

        return tube
