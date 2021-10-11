import numpy as np
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import rdmolops
from collections import deque
import itertools
from collections import deque
from rdkit.Chem import rdChemReactions
import networkx as nx
from rdkit.Chem import QED
from sascorer import *



def get_clogp_score(smiles, logp_m, logp_s, sascore_m, sascore_s, cycle_m, cycle_s):
	
	logp_value = Descriptors.MolLogP(MolFromSmiles(smiles))
	sascore = -sascorer.calculateScore(MolFromSmiles(smiles)) 
	cycle_score = 0


	cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(MolFromSmiles(smiles))))
	if len(cycle_list)==0:
		cycle_length = 0
	else:
		cycle_length = max([len(j) for j in cycle_list])
		if cycle_length <= 6:
			cycle_length = 0
		else:
			cycle_length = cycle_length - 6
	cycle_score = -cycle_length

	logP_value_normalized = (logp_value - logp_m)/logp_s
	sascore_normalized = (sascore - sascore_m)/sascore_s
	cycle_score_normalized = (cycle_score- cycle_m)/cycle_s

	return logP_value_normalized + sascore_normalized + cycle_score_normalized


def get_qed_score(smiles):
	score = QED.qed(rdkit.Chem.MolFromSmiles(smiles))
	return score


def check(rxn):
	molecule_nodes = rxn.molecule_nodes
	template_nodes = rxn.template_nodes
	root = template_nodes[0]
	queue = deque([root])
	root.depth = 0
	order ={}
	visited = set([root.id])
	node2smiles={}
	order[0]=[root.id]

	while len(queue) > 0:
		x = queue.popleft()
		for y in x.children:
			if len(y.children) == 0:
				continue
			template = y.children[0]
			if template.id not in visited:
				queue.append(template)
				visited.add(template.id)
				template.depth = x.depth + 1
				if template.depth not in order:
					order[template.depth] = [template.id]
				else:
					order[template.depth].append(template.id)
	#print(order)
	maxdepth = len(order) - 1
	for t in range(maxdepth, -1, -1):
		for template_id in order[t]:
			template_node = template_nodes[template_id]
			template = template_node.template
			
			reactants =[]
			for reactant_node in template_node.children:
				if len(reactant_node.children) == 0:
					reactant = reactant_node.smiles
					reactants.append(reactant)
					node2smiles[reactant_node.id] = reactant
					#print(reactant)

				else:
					reactant = node2smiles[reactant_node.id]
					reactants.append(reactant)
			possible_templates = reverse_template(template)
			#print(reactants, possible_templates)
			reacts = [Chem.MolFromSmiles(reactant) for reactant in reactants]
			product = None
			for template_ in possible_templates:
				try:
					rn = rdChemReactions.ReactionFromSmarts(template_)
					products = rn.RunReactants(reacts)
					a,b = template_.split(">>")
					n1 = len(a.split("."))
					n2 = len(reacts)
					if n1 != n2:
						#print("dkm", template_, reactants, product)
						#exit(1)
						return False

					if len(products) > 0:
						product = products[0][0]
						
						break
				except:
					return False
			if product == None:
				return False
			else:
				product_id = template_node.parents[0].id
				node2smiles[product_id] = Chem.MolToSmiles(product)
	return True


def filter_dataset(rxn_trees):
	return_rxns=[]
	for i, rxn in enumerate(rxn_trees):
		if check(rxn):
			#print(i, "OK")
		#else:
		#	print(i, "not OK")
			return_rxns.append(rxn)
	return return_rxns



def get_mol_from_smiles(smiles):
	mol = Chem.MolFromSmiles(smiles)
	if mol is None:
		return None
	Chem.Kekulize(mol)
	return mol

def get_smiles_from_mol(mol):
	return Chem.MolToSmiles(mol, kekuleSmiles=True)

def read_multistep_rxns(filename):
	synthetic_routes = []
	scores = []
	with open(filename, 'r') as reader:
		lines = reader.readlines()
		for i,line in enumerate(lines):

			full_rxn = []
			reactions = line.strip().split(' ')
			for reaction in reactions[:-1]:
				product, reactants, template = reaction.split('$')
				full_rxn.append([product, reactants, template])
				n1 = len(reactants.split("."))
				p1,p2 = template.split(">>")
				#n2 = len(p1.split("."))
				#if n1 !=n2:
				#	print(i, n1,n2, template, reaction)
			synthetic_routes.append(full_rxn)
			scores.append(float(reactions[-1]))
	return synthetic_routes, scores


def reverse_template(template):# for uspto
	p2, p1 = template.split(">>")
	p2 = p2.split(".")
	p2_list =  list(itertools.permutations(p2))
	reactant_list = [".".join(p2) for p2 in p2_list]

	return [">>".join([p2, p1])  for p2 in reactant_list]



def get_template_order(rxn):
	mol_nodes = rxn.molecule_nodes
	tem_nodes = rxn.template_nodes

	order={}
	root = tem_nodes[0]
	queue = deque([root])
	visisted = set([root.id])
	root.depth = 0
	order[0] =[root.id]
	
	while len(queue) > 0:
		x = queue.popleft()
		for y in x.children:
			if len(y.children) == 0: # starting molecule
				continue
			template = y.children[0] 
			if template.id not in visisted:
				queue.append(template)
				visisted.add(template.id)
				template.depth = x.depth + 1
				if template.depth not in order:
					order[template.depth] = [template.id]
				else:
					order[template.depth].append(template.id)
	return order


