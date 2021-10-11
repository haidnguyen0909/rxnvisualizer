
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Draw, AllChem

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as patches
import numpy as np
import networkx as nx
import os
from reaction import ReactionTree, extract_starting_reactants, StartingReactants, Templates, extract_templates,stats
from jinja2 import Template
from PIL import Image, ImageDraw
import tempfile
import subprocess


IMAGE_FOLDER = tempfile.mkdtemp()



def readFile(input_fn):
	synthetic_routes =[]
	scores=[]
	for raw_line in open(input_fn, "r").readlines():
		line = raw_line.strip()
		eles = line.strip().split(' ')
		score = -float(eles[-1])
		reactions = " ".join(eles[1:-1])
		synthetic_routes.append(reactions)
		scores.append(score)
	return synthetic_routes, scores



def crop_image(img, margin=20):
    x0_lim = img.width
    y0_lim = img.height
    x1_lim = 0
    y1_lim = 0
    for x in range(0, img.width):
        for y in range(0, img.height):
            if img.getpixel((x, y)) != (255, 255, 255):
                if x < x0_lim:
                    x0_lim = x
                if x > x1_lim:
                    x1_lim = x
                if y < y0_lim:
                    y0_lim = y
                if y > y1_lim:
                    y1_lim = y
    x0_lim = max(x0_lim, 0)
    y0_lim = max(y0_lim, 0)
    x1_lim = min(x1_lim + 1, img.width)
    y1_lim = min(y1_lim + 1, img.height)
    # Then crop to this area
    cropped = img.crop((x0_lim, y0_lim, x1_lim, y1_lim))
    # Then create a new image with the desired padding
    out = Image.new(
        img.mode,
        (cropped.width + 2 * margin, cropped.height + 2 * margin),
        color="white",
    )
    out.paste(cropped, (margin + 1, margin + 1))
    return out


def draw_rounded_rectangle(img, color, arc_size= 20):
    x0, y0, x1, y1 = img.getbbox()
    x1 -= 1
    y1 -= 1
    copy = img.copy()
    draw = ImageDraw.Draw(copy)
    arc_size_half = arc_size // 2
    draw.arc((x0, y0, arc_size, arc_size), start=180, end=270, fill=color)
    draw.arc((x1 - arc_size, y0, x1, arc_size), start=270, end=0, fill=color)
    draw.arc((x1 - arc_size, y1 - arc_size, x1, y1), start=0, end=90, fill=color)
    draw.arc((x0, y1 - arc_size, arc_size, y1), start=90, end=180, fill=color)
    draw.line((x0 + arc_size_half, y0, x1 - arc_size_half, y0), fill=color)
    draw.line((x1, arc_size_half, x1, y1 - arc_size_half), fill=color)
    draw.line((arc_size_half, y1, x1 - arc_size_half, y1), fill=color)
    draw.line((x0, arc_size_half, x0, y1 - arc_size_half), fill=color)
    return copy


def molecule_to_image(molecule, frame_color):
	mol = Chem.MolFromSmiles(molecule.smiles)
	img = Draw.MolToImage(mol)
	cropped_img = crop_image(img)
	final_img = draw_rounded_rectangle(cropped_img, frame_color)
	#final_img.save("dkm.png")
	return final_img


def reaction_to_image(template_node, frame_color):
	template = template_node.template
	rxn = AllChem.ReactionFromSmarts(template)
	img = Draw.ReactionToImage(rxn)
	cropped_img = crop_image(img)
	final_img = draw_rounded_rectangle(cropped_img, frame_color)
	#final_img.save("dkm.png")
	return final_img


def molecules_to_images(molecules, frame_colors):
	images=[]
	for idx, molecule in enumerate(molecules):
		img = molecule_to_image(molecule, frame_colors[idx])
		#img.save(str(idx)+".png")
		images.append(img)
	return images

def templates_to_images(template_nodes, frame_colors):
	images=[]
	for idx, template_node in enumerate(template_nodes):
		img = reaction_to_image(template_node, frame_colors[idx])
		#img.save(str(idx)+".png")
		images.append(img)
	return images

def save_molecule_images(molecule_nodes, frame_colors):
	global IMAGE_FOLDER
	images = molecules_to_images(molecule_nodes, frame_colors)
	spec ={}
	for idx in range(len(molecule_nodes)):
		molecule_node = molecule_nodes[idx]
		if molecule_node.id !=idx:
			print("IDX not consistent")
			exit(1)

		frame_color = frame_colors[idx]
		image_filepath = os.path.join(IMAGE_FOLDER, "mol_" + str(idx) +".png")
		image_obj = images[idx]
		image_obj.save(image_filepath)
		spec["mol_" + str(idx)] = image_filepath
	return spec

def save_template_images(template_nodes, frame_colors):
	global IMAGE_FOLDER
	images = templates_to_images(template_nodes, frame_colors)
	spec ={}
	for idx in range(len(template_nodes)):
		template_node = template_nodes[idx]
		frame_color = frame_colors[idx]
		image_filepath = os.path.join(IMAGE_FOLDER, "tem_" + str(idx) +".png")
		image_obj = images[idx]
		image_obj.save(image_filepath)
		spec["tem_" + str(idx)] = image_filepath
	return spec

def make_graphviz_image(mol_specs, tem_specs, edges):
	mol_frame_color ="green"
	tem_frame_color ="red"
	
	def _create_image():
		txt = template.render(molecules = mol_specs, templates = tem_specs, edges = edges)
		_, input_name = tempfile.mkstemp(suffix=".dot")
		_, output_name = tempfile.mkstemp(suffix=".png")
		with open(input_name, "w") as this_fileobj:
			this_fileobj.write(txt)
		cmd = "dot -Tpng " + input_name + " -o " + output_name
		subprocess.call(cmd.split(" "))
		return output_name


	template_filepath = os.path.join("", "reaction_template.dot")
	with open(template_filepath, "r") as fileobj:
		template = Template(fileobj.read())

	output_img = _create_image()
	return Image.open(output_img)


def route2image(route, score):

	reactions = route.strip().split(" ")
	reactions_array=[]
	for reaction in reactions:
		prod_reactant_template = reaction.split("$")
		product = prod_reactant_template[0]
		reactants = prod_reactant_template[1]
		template = prod_reactant_template[2]
		reactions_array.append([product, reactants, template])
	rxn = ReactionTree(reactions_array)
	tem_frame_color ="red"
	mol_frame_color ="green"
	tem_specs = save_template_images(rxn.template_nodes, [tem_frame_color]*len(rxn.template_nodes))
	mol_specs = save_molecule_images(rxn.molecule_nodes, [mol_frame_color] * len(rxn.molecule_nodes))

	prod_img = Image.open(mol_specs['mol_0'])
	startings = [Chem.MolFromSmiles(mol.smiles) for mol in rxn.molecule_nodes if len(mol.children)==0]

	edges =[]
	for id1, template_node in enumerate(rxn.template_nodes):
		parent = template_node.parents[0]
		id2 = parent.id
		edges.append(("mol_" + str(id2), "tem_" + str(id1)))
		for child in template_node.children:
			id2 = child.id
			edges.append(("tem_" + str(id1), "mol_" + str(id2)))

	rxn_img = make_graphviz_image(mol_specs, tem_specs, edges)
	return rxn_img, Draw.MolsToGridImage(startings, molsPerRow=10), prod_img, tem_specs






'''
def route2image(route, score):
	reactions = route.strip().split(" ")
	reactions_array=[]
	for reaction in reactions:
		prod_reactant_template = reaction.split("$")
		product = prod_reactant_template[0]
		reactants = prod_reactant_template[1]
		template = prod_reactant_template[2]
		reactions_array.append([product, reactants, template])
	rxn = ReactionTree(reactions_array)
	tem_frame_color ="yellow"
	tem_specs = save_template_images(rxn.template_nodes, [tem_frame_color]*len(rxn.template_nodes))
	#name2imgs ={}
	#for name, img in tem_specs.items():
	#	name2imgs[name]=Image.open(img)



	startings=[]
	for molecule_node in rxn.molecule_nodes:
		if len(molecule_node.children)==0:
			startings.append(molecule_node)
	

	
	startings = [Chem.MolFromSmiles(mol.smiles) for mol in startings]
	product = Chem.MolFromSmiles(rxn.molecule_nodes[0].smiles)

	molecule_nodes = rxn.molecule_nodes
	template_nodes = rxn.template_nodes
	root = molecule_nodes[0]
	output_img = make_graphviz_image(molecule_nodes, template_nodes)
	return output_img, Draw.MolsToGridImage(startings, molsPerRow=5), Draw.MolsToGridImage([product], molsPerRow=5, legends=[str(score)]), tem_specs
'''	
#route="O=S(=O)(Nc1ccc(F)cc1F)c1cnccc1Cl$O=S(=O)(Cl)c1cnccc1Cl.Nc1ccc(F)cc1F$([#17]-[S;H0;+0:2](=[*:1])(=[*:3])-[*:4]).([*:13]-[NH2;+0:12])>>([*:13]-[NH;+0:12]-[S;H0;+0:2](=[*:1])(=[*:3])-[*:4]) O=S(=O)(Cl)c1cnccc1Cl$O=S(=O)(O)c1cnccc1Cl.ClP(Cl)(Cl)(Cl)Cl$([#8]-[S;H0;+0:2](=[*:1])(=[*:3])-[*:4]).([#15]-[Cl;H0;+0:11])>>([*:1]=[S;H0;+0:2](=[*:3])(-[*:4])-[Cl;H0;+0:11]) Nc1ccc(F)cc1F$O=[N+]([O-])c1ccc(F)cc1F$([#8]-[N+;H0:1](=[#8])-[*:2])>>([*:2]-[NH2;+0:1]) O=[N+]([O-])c1ccc(F)cc1F$O=C(O)c1c(F)ccc([N+](=O)[O-])c1F$([*:7]:[c;H0;+0:8](:[*:9])-[C])>>([*:7]:[cH;+0:8]:[*:9])"

#route = "COCC(=O)Nc1nc2ccc(Br)cc2s1$COCC(=O)O.Nc1nc2ccc(Br)cc2s1$([#8]-[C;H0;+0:12](-[*:11])=[*:13]).([*:15]-[NH2;+0:14])>>([*:11]-[C;H0;+0:12](=[*:13])-[NH;+0:14]-[*:15]) COCC(=O)O$CO.O=C(O)CCl$([*:1]-[OH;+0:2]).([#17]-[CH2;+0:6]-[*:4])>>([*:1]-[O;H0;+0:2]-[CH2;+0:6]-[*:4]) Nc1nc2ccc(Br)cc2s1$NC(=S)Nc1ccc(Br)cc1$([*:1]-[C;H0;+0:2](=[S;H0;+0:3])-[NH;+0:4]-[*:5]:[cH;+0:6]:[*:7])>>([*:1]-[c;H0;+0:2]1:[n;H0;+0:4]:[*:5]:[c;H0;+0:6](:[*:7]):[s;H0;+0:3]:1) CO$COc1cc2nc(Cl)nc(NN)c2cc1OC$([*:1]-[O;H0;+0:2]-[c])>>([*:1]-[OH;+0:2]) NC(=S)Nc1ccc(Br)cc1$Nc1nc2ccc(Br)cc2s1$([*:1]-[c;H0;+0:2]1:[n;H0;+0:3]:[*:4]:[c;H0;+0:9](:[*:8]):[s;H0;+0:10]:1)>>([*:1]-[C;H0;+0:2](=[S;H0;+0:10])-[NH;+0:3]-[*:4]:[cH;+0:9]:[*:8]"
#route2image(route, 0.66)
#img.save("best3.png")
#for starting in startings:
#	starting.save("dkm1.png")










