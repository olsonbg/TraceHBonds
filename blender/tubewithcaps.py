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
                                                  radius=self.radius)
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

    def merge_verts(self, tube):
        bm = bmesh.new()

        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(tube.data)  # Fill bmesh with object mesh

        # get vertices of merged tube and caps which have 3 connections. They
        # should be on the ends of the tube, and the caps
        tube_end_verts = self.verts_withNconnections(bm, 3)

        # TODO: Do not hardcode number (0.001). Find a way to figure out
        # a good value for this.
        bmesh.ops.remove_doubles(bm, verts=tube_end_verts, dist=0.001)

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

        #  for v in verts:
        #      print(v.co)

        # Get the coordinates of the the vertices
        verts_co = [(ob.matrix_world @ v.co) for v in verts]

        return verts_co

    def vertsCenter(self, verts):
        # Get geometric center of cap and tube ends.
        vSum = mathutils.Vector((0, 0, 0))
        for v in verts:
            vSum += v

        return vSum/len(verts)

    def extrude_tube(self):

        polyline, curvedata = self.make_polyline()

        bevel_ob = self.make_bevelobject()
        curvedata.bevel_mode = 'OBJECT'
        curvedata.bevel_object = bevel_ob


        # Convert to a mesh
        obj = bpy.context.scene.objects.get(self.name)
        obj.select_set(True)

        selection_names = [obj.name for obj in bpy.context.selected_objects]

        #  bpy.context.scene.objects[self.name].select_set(True)
        #  bpy.data.objects[self.name].select = True
        bpy.ops.object.convert(target='MESH')

        # Done with bevel object, so delete it.
        bpy.ops.object.select_all(action='DESELECT')
        bevel_ob.select_set(True)
        bpy.ops.object.delete()

        return polyline
        #  return bpy.data.objects[self.name]

    def make_polyline(self):

        curvedata = bpy.data.curves.new(name=self.pathName, type='CURVE')
        curvedata.dimensions = '3D'

        objectdata = bpy.data.objects.new(self.name, curvedata)
        objectdata.location = (0, 0, 0)  # object origin
        # API Change: New objects must be linked to a collection, not scene.
        bpy.context.view_layer.active_layer_collection.collection.objects.link(objectdata)

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

        cap_end_verts_co = self.verts_withNconnections_co(cap, 3)
        cn = mathutils.geometry.normal(*cap_end_verts_co)

        if cap_normal == mathutils.Vector((0, 0, 1)):
            cn = -1.0 * cn


        # Get normal to tube end. Pick three points equally spaced to determine
        # normals
        one_third = int(len(tube_end_co)/3)
        end = [tube_end_co[0],
               tube_end_co[one_third],
               tube_end_co[2*one_third]]

        normal = mathutils.geometry.normal(*end)
        # normal_sphere = mathutils.Vector((0, 0, -1))
        to_rotate = cap_normal.rotation_difference(normal)
        to_rotate = cn.rotation_difference(normal)
        cap.rotation_mode = 'QUATERNION'
        # Rotate the cap so that it is at the same angle as the tube ends
        cap.rotation_quaternion = to_rotate

        # Apply transforms for later calculations.
        bpy.ops.object.transform_apply(scale=True)
        


        # Get new cap verts
        capverts = self.verts_withNconnections_co(cap, 3)

        capCenter = self.vertsCenter(capverts)
        tubeCenter = self.vertsCenter(tube_end_co)

        # Rotate cap to align vertices with tube vertices. The vertices will
        # be very close to aligned, but not perfectly.
        v1 = capverts[0]    - capCenter
        v2 = tube_end_co[0] - tubeCenter

        #  if cap_normal == mathutils.Vector((0, 0, 1)):
        #      v1 = -1.0 * v1

        anglediff = v1.rotation_difference(v2)

        cap.delta_rotation_quaternion = -1.0*anglediff

        # Move cap to the tube end
        cap.location = location

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

        # Make hemispheres with same number of segments as the tube ends
        capStart = self.makeHemisphere(self.startCapName, num_verts)
        capEnd = self.makeHemisphere(self.endCapName, num_verts)

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
        capEnd.select_set(True)
        tube.select_set(True)
        capStart.select_set(True)
        bpy.context.view_layer.objects.active = tube
        #  bpy.context.scene.objects.active = tube

        bpy.ops.object.join()

        self.merge_verts(tube)

        return tube
