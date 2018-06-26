#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Name:        hbhcain.py
# Purpose:     Show hydrogen bonded chains in blender using the output from
#              TraceHBonds --json. An animation is also setup, showing a
#              rotating PBC box. Uses tubewithcaps.py for generating tubes, and
#              diverging_map.py for coloring tubes by length with a divergent
#              color map.
#
# Author:      Brian G. Olson (olsonbg@gmail.com)
#
# Created:     27 October 2015
# Copyright:   (c) 2015
#
# Usage:       blender -P hbchain.py -- -f HBonds1.json
#
# For help:    blender -P hbchain.py -- --help
# -----------------------------------------------------------------------------
#

import numpy as np
import json
import bpy
import mathutils
import math
# Import from the same directory this script is run from
import os
import sys
sys.path.append(os.path.dirname(__file__))

from tubewithcaps import TubeWithCaps
from diverging_map import ColorMapCreator

import time
start_time = time.time()

NumFrames = 120
Scale = 5.0
Radius = 0.1
WhiteIndex = 3


#
# used for testing
#
def put_marker(co, name="Test"):
    obj_empty = bpy.data.objects.new(name, None)
    bpy.context.scene.objects.link(obj_empty)
    obj_empty.location = co


# Add PBC box, camera, and light source
def draw_box(box, scale):
    PBC = box["xyz"]
    PBC_2 = [x/2.0/scale for x in PBC]

    # Empty object used for rotation/animation, placed at center of PBC box,
    # and will also be the parent of the camera and lamps.
    bpy.ops.object.empty_add(type='CUBE', radius=0.1, location=(PBC_2))
    bpy.context.object.name = 'RotationBox'
    rotation_box = bpy.data.objects['RotationBox']

    bpy.ops.mesh.primitive_cube_add(radius=PBC_2[0], location=(PBC_2))
    bpy.ops.object.modifier_add(type='WIREFRAME')
    bpy.context.object.modifiers["Wireframe"].use_relative_offset = True
    bpy.context.object.modifiers["Wireframe"].use_even_offset = False
    bpy.context.object.name = 'PBC'

    # Add camera, looking down the z axis. RotationBox will be the parent, so
    # location is relative to 'RotationBox'
    camera_co = (0.0, 2.0*PBC_2[1], 10.0*PBC_2[2])
    # Angle the camera to it looks at the center of the PBC box (RotationBox)
    angle = -1.0 * math.atan(camera_co[1]/camera_co[2])
    bpy.ops.object.camera_add(view_align=True,
                              location=camera_co,
                              rotation=(angle, 0, 0))
    bpy.context.object.name = 'Camera'
    bpy.data.objects['Camera'].parent = rotation_box

    # Add light source

    # Create new lamp datablock
    lamp_data = bpy.data.lamps.new(name="lamp-center", type='SUN')

    # Create new object with our lamp datablock
    lamp_object = bpy.data.objects.new(name="lamp-center",
                                       object_data=lamp_data)
    lamp_data.use_nodes = True

    #  Link lamp object to the scene so it'll appear in this scene
    bpy.context.scene.objects.link(lamp_object)

    lamp_object.location = (0, 0, 10.1 * PBC_2[2])

    bpy.data.objects['lamp-center'].parent = rotation_box

    nodes = bpy.data.lamps['lamp-center'].node_tree.nodes
    # clear all nodes to start clean
    for node in nodes:
        nodes.remove(node)

    # create emission node
    node_emission = nodes.new(type='ShaderNodeEmission')
    node_emission.inputs[0].default_value = (1, 1, 1, 1)  # While RGBA
    node_emission.inputs[1].default_value = 3.0  # Strength
    node_emission.location = 0, 0

    # create output node
    node_output = nodes.new(type='ShaderNodeOutputMaterial')
    node_output.location = 200, 0

    # Link nodes
    links = bpy.data.lamps['lamp-center'].node_tree.links
    links.new(node_emission.outputs[0], node_output.inputs[0])


def element_sphere(obname, element, scale):
    # Typical CPK coloring assignments from
    # https://en.wikipedia.org/wiki/CPK_coloring
    element_colors = {
                      'H' : (221, 221, 221, 255),
                      'C' : ( 34,  34,  34, 255),
                      'N' : (135, 206, 235, 255),
                      'O' : (255,  13,  13, 255),
                      'F' : ( 31, 240,  31, 255),
                      'Cl': ( 31, 240,  31, 255),
                      'Br': (153,  34,   0, 255),
                      'I' : (102,   0, 187, 255),
                      'He': (  0, 255, 255, 255),
                      'Ne': (  0, 255, 255, 255),
                      'Ar': (  0, 255, 255, 255),
                      'Xe': (  0, 255, 255, 255),
                      'Kr': (  0, 255, 255, 255),
                      'P' : (255, 153,   0, 255),
                      'S' : (221, 221,   0, 255),
                      'other' : (221, 119, 255, 255)
                     }  # noqa

    # Van der Waals Radius of the elements (angstroms)
    element_sizes = {
                     'H' : 1.20,
                     'Zn': 1.39,
                     'He': 1.39,
                     'Cu': 1.40,
                     'F' : 1.47,
                     'O' : 1.52,
                     'Ne': 1.54,
                     'N' : 1.55,
                     'Au': 1.66,
                     'C' : 1.70,
                     'Ag': 1.72,
                     'Cl': 1.75,
                     'P' : 1.80,
                     'S' : 1.80,
                     'other' : 1.70
                    }  # noqa

    mat_name = 'ElementMaterial-{0:s}'.format(element)

    mat = None
    if mat_name not in bpy.data.materials:
        mat = bpy.data.materials.new(mat_name)
        if element in element_colors:
            color = element_colors[element]
        else:
            color = element_colors['other']
            print("No color definition for {0:s}.".format(element))

        color = [x/255.0 for x in color]

        material_nodes(mat, color)
    else:
        mat = bpy.data.materials[mat_name]

    # Sphere for this element
    bpy.ops.object.select_all(action='DESELECT')

    if element in element_sizes:
        e_size = element_sizes[element]
    else:
        e_size = element_sizes['other']
        print("No size definition for {0:s}.".format(element))

    # For ball and stick divide by 5.0
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0, 0, 0),
                                         size=e_size/scale/5.0)
    sphere = bpy.context.object
    sphere.name = obname
    bpy.ops.object.shade_smooth()
    sphere.data.materials.append(mat)
    sphere.select = False

    bpy.context.scene.objects.unlink(sphere)
    return(sphere)


# Cylinder for bonds
def element_cylinder(name, bondradius, scale):
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.mesh.primitive_cylinder_add(location=(0, 0, 0),
                                        radius=bondradius/scale)
    cylinder = bpy.context.object
    cylinder.name = name
    cylinder.rotation_mode = 'QUATERNION'
    cylinder.data.materials.append(bpy.data.materials['BondMaterial'])
    cylinder.select = False

    bpy.context.scene.objects.unlink(cylinder)

    return(cylinder)


def draw_structure(data, scale):
    # Store objects here, then link them to the scene at the end of this
    # function.
    obs = list()

    # Empty object used as parent for chemical structure.
    ob_structure = bpy.data.objects.new('Structure', None)
    ob_structure.hide        = False
    ob_structure.hide_render = False
    ob_structure.hide_select = False
    obs.append(ob_structure)

    # Materials
    mat = bpy.data.materials.new("BondMaterial")
    material_nodes(mat, (0.02, 0.02, 0.02, 1))

    # Initial cylinder for bonds.
    cylinder_axis = mathutils.Vector((0, 0, 1))
    bondradius = 1.2/5.0/4.0
    v_double = 1.5*bondradius/scale
    v_triple = 3.0*bondradius/scale

    ob_elements = {}
    ob_bondtypes = {}
    Ats = []

    for mol in data:
        ob_molecule = None
        ob_atoms = None
        ob_element = None
        ob_bonds = None
        ob_bondtype = None

        if 'name' in mol:
            # Empty object used as parent for molecules.
            ob_molecule = bpy.data.objects.new(mol['name'], None)
            ob_molecule.parent = ob_structure
            obs.append(ob_molecule)
            print("Molecule: {0:s}".format(mol['name']))

        if 'atoms' in mol:
            start_time2 = time.time()
            # Empty object used as parent for atoms.
            atom_parent_name = 'Atoms ({0:s})'.format(mol['name'])
            ob_atoms = bpy.data.objects.new(atom_parent_name, None)
            ob_atoms.parent = ob_molecule
            obs.append(ob_atoms)

            for atom in mol['atoms']:
                loc = [x/scale for x in atom['location']]

                # element = "{0:s} {1:s}".format(atom['element'],
                #                                atom['forcefield'])
                atom_name = "{0:s} {1:s} ({2:s})".format(atom['element'],
                                                         atom['forcefield'],
                                                         mol['name'])

                if atom_name not in ob_elements:
                    sphere_name = 'Sphere-{0:s}'.format(atom_name)
                    ob_elements[atom_name] = element_sphere(sphere_name,
                                                            atom['element'],
                                                            scale)

                    # Empty element to be parent of these atoms
                    ob_element = bpy.data.objects.new(atom_name, None)
                    ob_element.parent = ob_atoms
                    obs.append(ob_element)
                else:
                    ob_element = bpy.data.objects[atom_name]

                ob = ob_elements[atom_name].copy()
                ob.location = loc
                ob.parent = ob_element
                ob.name = atom['name']

                obs.append(ob)

                # Used for determining bond locations
                Ats += [None]*len(mol['atoms'])
                Ats[atom['ID']] = [loc, atom['forcefield'], atom['name']]

            print("Atoms: %s seconds" % (time.time() - start_time2))

        if 'bonds' in mol:
            start_time2 = time.time()
            # Empty object used as parent for bonds
            bond_parent_name = 'Bonds ({0:s})'.format(mol['name'])
            ob_bonds = bpy.data.objects.new(bond_parent_name, None)
            ob_bonds.parent = ob_molecule
            obs.append(ob_bonds)

            for bond in mol['bonds']:
                ID1 = bond["IDs"][0]
                ID2 = bond["IDs"][1]

                v1 = mathutils.Vector(Ats[ID1][0])
                v2 = mathutils.Vector(Ats[ID2][0])
                vbond = v2-v1
                height = vbond.magnitude
                loc = v1+vbond/2.0
                to_rotate = cylinder_axis.rotation_difference(vbond)

                bond_name = "{0:s}-{1:s} [{2:d}] ({3:s})".format(Ats[ID1][1],
                                                                 Ats[ID2][1],
                                                                 bond['order'],
                                                                 mol['name'])
                if bond_name not in ob_bondtypes:
                    cylinder_name = 'Cylinder-{0:s}'.format(bond_name)
                    ob_bondtypes[bond_name] = element_cylinder(cylinder_name,
                                                               bondradius,
                                                               scale)

                    # Empty element to be parent of these bonds
                    ob_bondtype = bpy.data.objects.new(bond_name, None)
                    ob_bondtype.parent = ob_bonds
                    obs.append(ob_bondtype)
                else:
                    ob_bondtype = bpy.data.objects[bond_name]

                ob = ob_bondtypes[bond_name].copy()
                ob.parent = ob_bondtype
                ob.dimensions[2] = height
                ob.rotation_quaternion = to_rotate

                if bond['order'] == 1:
                    ob.location = loc
                    ob.name = '{0:d}-{1:d}'.format(ID1, ID2)
                    obs.append(ob)
                elif bond['order'] == 2:
                    # any unit vector perpendicular to vbond (the bond vector)
                    # will suffice.
                    orthogonal = vbond.orthogonal().normalized()

                    ob.location = loc + v_double * orthogonal

                    ob2 = ob.copy()
                    ob2.location = loc - v_double * orthogonal

                    ob.name  = '{0:d}-{1:d} b1'.format(ID1, ID2)
                    ob2.name = '{0:d}-{1:d} b2'.format(ID1, ID2)

                    obs.append(ob)
                    obs.append(ob2)
                elif bond['order'] == 3:
                    # any unit vector perpendicular to vbond (the bond vector)
                    # will suffice.
                    orthogonal = vbond.orthogonal().normalized()

                    ob.location = loc

                    ob2 = ob.copy()
                    ob2.location = ob.location + v_triple * orthogonal

                    ob3 = ob.copy()
                    ob3.location = ob.location - v_triple * orthogonal

                    ob.name  = '{0:d}-{1:d} b1'.format(ID1, ID2)
                    ob2.name = '{0:d}-{1:d} b2'.format(ID1, ID2)
                    ob3.name = '{0:d}-{1:d} b3'.format(ID1, ID2)

                    obs.append(ob)
                    obs.append(ob2)
                    obs.append(ob3)

            print("Bonds: %s seconds " % (time.time() - start_time2))

    # Cleanup
    bpy.ops.object.select_all(action='DESELECT')
    #  Select objects to delete
    for ob_name in ob_bondtypes:
        ob_bondtypes[ob_name].select = True

    for ob_name in ob_elements:
        ob_elements[ob_name].select = True

    bpy.ops.object.delete()

    print("Linking objects to scene")
    for ob in obs:
        bpy.context.scene.objects.link(ob)
    print("Updating scene")

    bpy.context.scene.update()


def draw_hbchain(chain, i, radius, scale):
    coordinates = list()
    numatoms = len(chain)

    for atom in chain:
        coordinates.append([x/scale for x in atom["location"]])

    tubeName   = 'Chain{0:04d}'.format(i)
    groupName  = 'GroupLength{0:02d}'.format(numatoms)
    lengthName = 'Length{0:02d}'.format(numatoms)

    # If groupName doesn't exist, then neither does the empty LengthName
    if groupName not in bpy.data.groups:
        bpy.ops.group.create(name=groupName)
        # Make empty as parent for chains of this length
        bpy.ops.object.empty_add(type='CUBE', radius=0.1, location=(0, 0, 0))
        bpy.context.object.name = lengthName
        bpy.context.object.parent = bpy.data.objects['Chains']

    group     = bpy.data.groups[groupName]
    ob_length = bpy.data.objects[lengthName]

    tube = TubeWithCaps(tubeName, radius, coordinates)
    ob_tube = tube.make()

    ob_tube.parent = ob_length
    group.objects.link(ob_tube)

    bpy.ops.object.select_all(action='DESELECT')
    return numatoms


def material_nodes(mat, color):
    mat.use_nodes = True
    nodes = mat.node_tree.nodes

    # clear all nodes to start clean
    for node in nodes:
        nodes.remove(node)

    node_rgba = nodes.new(type='ShaderNodeRGB')
    node_rgba.outputs[0].default_value = color
    node_rgba.location = (-200, 400)

    node_rgbaBack = nodes.new(type='ShaderNodeRGB')
    node_rgbaBack.outputs[0].default_value = (0, 0, 0, 1)
    node_rgbaBack.location = (-200, 200)

    node_mix = nodes.new(type='ShaderNodeMixRGB')
    node_mix.inputs[0].default_value = 0.7
    node_mix.location = (0, 200)

    node_layerweight = nodes.new(type='ShaderNodeLayerWeight')
    node_layerweight.inputs[0].default_value = 0.5
    node_layerweight.location = (200, 500)

    node_bsdf1 = nodes.new(type='ShaderNodeBsdfDiffuse')
    node_bsdf1.inputs[1].default_value = 0.0
    node_bsdf1.location = (200, 350)

    node_bsdf2 = nodes.new(type='ShaderNodeBsdfDiffuse')
    node_bsdf2.inputs[1].default_value = 0.0
    node_bsdf2.location = (200, 200)

    node_fresnel = nodes.new(type='ShaderNodeFresnel')
    node_fresnel.inputs[0].default_value = 1.450
    node_fresnel.location = (400, 500)

    node_mixshader1 = nodes.new(type='ShaderNodeMixShader')
    node_mixshader1.location = (400, 350)

    node_glossy = nodes.new(type='ShaderNodeBsdfGlossy')
    node_glossy.inputs[0].default_value = (1, 1, 1, 1)
    node_glossy.inputs[1].default_value = 0.050
    node_glossy.location = (400, 200)

    node_mixshader2 = nodes.new(type='ShaderNodeMixShader')
    node_mixshader2.location = (600, 350)

    node_output = nodes.new(type='ShaderNodeOutputMaterial')
    node_output.location = (800, 350)

    links = mat.node_tree.links
    links.new(node_rgba.outputs[0], node_mix.inputs[1])
    links.new(node_rgba.outputs[0], node_bsdf1.inputs[0])
    links.new(node_rgbaBack.outputs[0], node_mix.inputs[2])
    links.new(node_mix.outputs[0], node_bsdf2.inputs[0])
    links.new(node_layerweight.outputs[1], node_mixshader1.inputs[0])
    links.new(node_bsdf1.outputs[0], node_mixshader1.inputs[1])
    links.new(node_bsdf2.outputs[0], node_mixshader1.inputs[2])
    links.new(node_fresnel.outputs[0], node_mixshader2.inputs[0])
    links.new(node_mixshader1.outputs[0], node_mixshader2.inputs[1])
    links.new(node_glossy.outputs[0], node_mixshader2.inputs[2])
    links.new(node_mixshader2.outputs[0], node_output.inputs[0])


def generate_colormap(num_low, num_high):
    # ColorMap for chain lengths
    RGB1  = np.array([ 59,  76, 192])
    RGB2  = np.array([180,   4,  36])
    WHITE = np.array([221, 221, 221])

    colormap1 = ColorMapCreator(RGB1,  WHITE, numColors=num_low,    quiet=True)
    colormap2 = ColorMapCreator(WHITE, RGB2,  numColors=num_high+1, quiet=True)

    ColorMap = colormap1.colorMap.tolist()[:] + colormap2.colorMap.tolist()[1:]

    for i, c in enumerate(ColorMap):
        c.append(1)

    return ColorMap


# Make materials for objects:
#   1 per group
#   1 for PBC box
#
def generate_materials(white_index):
    import re

    # The number of materials to generatre depends on the maximum chainlength
    mat_count = 0
    for group in bpy.data.groups:
        if group.name.startswith("GroupLength"):
            regex = re.compile(r'\d+')
            j = int(regex.search(group.name).group(0))
            i = int((j - 1)/2)
            if i > mat_count:
                mat_count = i

    if mat_count < white_index:
        colormap = generate_colormap(mat_count, 0)
    else:
        colormap = generate_colormap(white_index, mat_count-white_index)

    # Generate textures
    for i in range(mat_count):
        mat_name = "GroupLength{0:02d}".format(2*i+3)
        mat = bpy.data.materials.new(mat_name)
        material_nodes(mat, colormap[i])

    # Also add a material for the PBC box
    mat = bpy.data.materials.new("PBCMaterial")
    mat.diffuse_color = 0.0, 0.0, 0.0
    # Don't want the PBC box to shine much at all when directly in the light
    mat.specular_intensity = 0.02


#
# Assign generated materials to objects
#
def assignMaterials():
    # Assign a texture to each chain length group

    for group in bpy.data.groups:
        mat = bpy.data.materials[group.name]
        #  print("Group: {0:s}".format(group.name))
        for object in group.objects:
            #  print("  Chain: {0:s}".format(object.name))
            if object.data is not None:
                if len(object.material_slots) < 1:
                    object.data.materials.append(mat)
                else:
                    object.material_slots[0].material = mat

    # Give the PBC box its material
    ob = bpy.data.objects["PBC"]
    mat = bpy.data.materials["PBCMaterial"]
    if len(ob.material_slots) < 1:
        ob.data.materials.append(mat)
    else:
        ob.material_slots[0].material = mat
    ob.modifiers["Wireframe"].material_offset = 0


#
# Animation settings
#
def setupAnimation(numframes):
    # bpy.ops.object.select_all(action='DESELECT')
    rotbox = bpy.data.objects['RotationBox']

    bpy.context.scene.frame_start = 1
    bpy.context.scene.frame_end   = numframes

    # Initial rotation
    rotbox.rotation_euler[0] = 0
    rotbox.rotation_euler[1] = 0
    rotbox.rotation_euler[2] = 0
    # Set this rotation as initial keyframe
    rotbox.keyframe_insert("rotation_euler", frame=1)

    # Final rotation
    rotbox.rotation_euler[0] = 0
    rotbox.rotation_euler[1] = 6.28319
    rotbox.rotation_euler[2] = 0
    # Set this rotation as keyframe
    rotbox.keyframe_insert("rotation_euler", frame=numframes + 1)

    # Use linear interpolation along frames
    for fcur in rotbox.animation_data.action.fcurves:
        fcur.extrapolation = 'LINEAR'


def main():
    import sys       # to get command line args
    import argparse  # to parse options for us and print a nice help message

    # get the args passed to blender after "--", all of which are ignored by
    # blender so scripts may receive their own arguments
    argv = sys.argv

    if "--" not in argv:
        argv = []  # as if no args are passed
    else:
        argv = argv[argv.index("--") + 1:]  # get all args after "--"

    # When --help or no args are given, print this help
    usage_text = \
        "Run blender with this script:\n\n" \
        "  blender --python " + __file__ + " -- [options]"

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=usage_text)

    # Possible types are: string, int, long, choice, float and complex.
    parser.add_argument("-f", "--file", dest="filename", type=str,
                        metavar='FILE', required=True,
                        help="This file will be loaded")

    parser.add_argument("-c", "--chemfile", dest="chemfile", type=str,
                        metavar='FILE', required=False, default=None,
                        help="This file describes the full chemical structure")

    parser.add_argument("-r", "--radius", dest="radius", type=float,
                        required=False, default=Radius,
                        help="Radius of the tubes. (default: %(default)f)")

    parser.add_argument("-n", "--numframes", dest="frames", type=int,
                        required=False, default=NumFrames,
                        help="Number of frames rendered for animation. (default: %(default)d)")

    parser.add_argument("-s", "--scale", dest="scale", type=float,
                        required=False, default=Scale,
                        help="Scale coordinates down by this factor. (default: %(default)f)")

    parser.add_argument("-w", "--white", dest="white_index", type=int,
                        required=False, default=WhiteIndex,
                        help="Location of 'white' in the divergent colormap. From the shortest chain to WHITE_INDEX the color will go blue to white. From WHITE_INDEX to the longest chain the color will go white to red. (default: %(default)d)")

    args = parser.parse_args(argv)

    if not argv:
        parser.print_help()
        return

    # Load data
    print("Loading data.")
    json_data = open(args.filename).read()
    hb_data = json.loads(json_data)

    scene = bpy.context.scene

    # Use Cycles render engine
    scene.render.engine = 'CYCLES'
    scene.cycles.samples = 200
    scene.cycles.preview_samples = 20
    scene.cycles.blur_glossy = 0.50
    scene.world.horizon_color = (101/255.0, 123/255.0, 131/255.0)

    # Start with an empty scene
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.object.select_all(action='TOGGLE')
    bpy.ops.object.delete(use_global=False)

    # Remove any materials which may be left
    for material in bpy.data.materials:
        bpy.data.materials.remove(material, do_unlink=True)

    # Empty object used as parent for chains.
    bpy.ops.object.empty_add(type='CUBE', radius=0.1, location=(0, 0, 0))
    bpy.context.object.name = 'Chains'

    print("Drawing hydrogen bond chains.")
    i = 0
    for data in hb_data:
        if "hbchain" in data:
            #  print("Drawing chain {0:d}.".format(i))
            draw_hbchain(data["hbchain"], i, args.radius, args.scale)
            i = i+1
        elif "PBC" in data:
            # Add PBC box, camera, and light source
            draw_box(data["PBC"], args.scale)

    print("Number of chains: {0:d}.".format(i))

    generate_materials(args.white_index)
    assignMaterials()

    if args.chemfile is not None:
        # Load chemical structure data
        print("Loading data.")
        json_data = open(args.chemfile).read()
        chem_data = json.loads(json_data)
        print("Drawing chemical structure.")
        for data in chem_data:
            if "molecules" in data:
                # Draw to full chemical structure
                draw_structure(data["molecules"], args.scale)

    setupAnimation(args.frames)

    # Turn off relationship lines in all views
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            area.spaces[0].show_relationship_lines = False

    bpy.data.objects["RotationBox"].hide        = True
    bpy.data.objects["RotationBox"].hide_render = True


if __name__ == "__main__":
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
