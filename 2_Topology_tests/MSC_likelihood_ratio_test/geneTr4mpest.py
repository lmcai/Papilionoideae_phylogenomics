from ete3 import Tree

x=open('G1456.trees').readlines()
out=open('G1456.4mpest.trees','a')


outgr=['Cercis_canadensis_ge','FC230_c98','FC232_c98']
outgr2=['TX04','HM2443_c98','RB56','G007_c98','HM2584_c98']
_exhausted  = object()

for l in x:
	t=Tree(l)
	sp=[node.name for node in t]
	outsp=[i for i in sp if i in outgr]
	if len(outsp)>0:
		outnode=t.get_monophyletic(values=outsp, target_attr="name")
		#check if outnode is empty
		if next(outnode, _exhausted) is _exhausted:
			t.set_outgroup(t&outsp[0])
		else:
			outnode=t.get_monophyletic(values=outsp, target_attr="name")
			for node in outnode:
				t.set_outgroup(node)
		d=out.write(t.write(format=9)+'\n')
	else:
		outsp=[i for i in sp if i in outgr2]
		outnode=t.get_monophyletic(values=outsp, target_attr="name")
		#check if outnode is empty
		if next(outnode, _exhausted) is _exhausted:
			t.set_outgroup(t&outsp[0])
		else:
			outnode=t.get_monophyletic(values=outsp, target_attr="name")
			for node in outnode:
				t.set_outgroup(node)
		d=out.write(t.write(format=9)+'\n')


out.close()
