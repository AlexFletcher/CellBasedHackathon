import os, sys
import cPickle as pickle
from sa_oa.tissueshape import centroid
from sa_oa.container import ordered_pids
from model_utils.celltissue_util import *

from model_utils.db_utilities import get_mesh,get_graph

from numpy import *
from pylab import *

import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection,PolyCollection
import matplotlib.pyplot as plt


from math import atan, cos,pi,sin,sqrt
from datetime import datetime
import time
a=datetime.now()

from sa_oa.ext.color import JetMap, GrayMap, GreenMap

def plot_pylab_figure(db,iteration,savename,prop_list,cmaps):
	mesh=get_mesh(db)
	graph=get_graph(db)
	db_position = db.get_property('position')
	poslist=list(db_position.itervalues())
	xvals,yvals = [[z[i] for z in poslist] for i in (0,1)]
	norm=max(max(xvals)-min(xvals),max(yvals)-min(yvals))
	position={}

	
	vindex = db.get_property('vindex')
	wall=db.get_property('wall')
	if mesh.degree()==3:
		zvals = [z[2] for z in poslist]
		for pid,pos in db_position.iteritems():
			position[pid]=((pos[0]-min(xvals))/norm,(pos[1]-min(yvals))/norm,pos[2])
		cpos={}
		for cid in mesh.wisps(3):
			cpos[cid]=[]
			for wid in mesh.borders(3,cid):
				for lid in mesh.borders(2,wid):
					cpair=[]
					for pid in mesh.borders(1,lid):
						cpair.append(position[int(pid)])
					cpos[cid].append((cpair[0],cpair[1]))
			cpos[cid]=list(set(cpos[cid]))
		ctc={}
		cell_to_coordinates={}
		
		octc={}
		for cid,pairs in cpos.iteritems():
			ctc[cid]=[]
			pv=[]
			for pair in pairs:
				pv.append(pair[0][2])
				pv.append(pair[1][2])
			for pair in pairs:
				if pair[0][2]==min(pv) and pair[1][2]==min(pv):
					ctc[cid].append(pair)
			octc[cid]=[]
			octc[cid].append(ctc[cid][0][0])
			octc[cid].append(ctc[cid][0][1])
			last=octc[cid][-1]
			while last!=octc[cid][0]:
				for pair in ctc[cid][1:]:
					if pair[0]==last:
						octc[cid].append(pair[1])
					elif pair[1]==last:
						octc[cid].append(pair[0])
					last=octc[cid][-1]
			octc[cid]=octc[cid][0:-1]
			cell_to_coordinates[cid]=[]
			for coord in octc[cid]:
				cell_to_coordinates[cid].append((coord[0],coord[1]))
		edge_coords={}
		wall_coords={}
		off=0.0001
		for eid in graph.edges():
			edge_coords[eid]=[]
			owid=wall[eid]			
			sid=graph.source(eid)
			tid=graph.target(eid)
			if vindex[sid]==vindex[tid]:
				if vindex[sid]==max(vindex.itervalues()):
					for wid in mesh.borders(3,sid):
						for cid in mesh.regions(2,wid):
							if vindex[cid]<vindex[sid]:
								centre=centroid(mesh,position,2,int(wid))
								basewall=int(wid)

				else:
					for wid in mesh.borders(3,sid):
						for cid in mesh.regions(2,wid):
							if vindex[cid]>vindex[sid]:
								centre=centroid(mesh,position,2,int(wid))
								basewall=int(wid)
						
				edgewall=wall[eid]
				lines=list(set(mesh.borders(2,edgewall))&set(mesh.borders(2,basewall)))
				line=int(lines[0])
				
				for pid in mesh.borders(1,line):
					pidpos=[position[int(pid)][0],position[int(pid)][1]]
					edge_coords[eid].append(pidpos)
					n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
					inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
					edge_coords[eid].append(inpos)
				edge_coords[eid][3]=edge_coords[eid][2]
				edge_coords[eid][2]=inpos

				if owid in list(wall_coords.iterkeys()):
					wall_coords[owid].append(edge_coords[eid][3])
					wall_coords[owid].append(edge_coords[eid][2])
					wall_coords[owid].append(edge_coords[eid][1])
				else:
					wall_coords[owid]=[edge_coords[eid][0],edge_coords[eid][1],edge_coords[eid][2]]
			else:
				centre=centroid(mesh,position,2,wall[eid])
				sid=graph.source(eid)
				wall_coords[owid]=[]
				for pidpos in cell_to_coordinates[sid]:
					n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
					inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
					edge_coords[eid].append(inpos)
					wall_coords[owid].append(inpos)
		
		for wid in mesh.wisps(2):
			if wid not in list(wall_coords.iterkeys()):
				wall_coords[wid]=[]
				centre=centroid(mesh,position,2,wid)
				if centre[2]==min(zvals):
					for pidpos in cell_to_coordinates[sid]:
						n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
						inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
						wall_coords[wid].append(inpos)
				elif centre[2]==max(zvals):
					for pidpos in cell_to_coordinates[sid]:
						n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
						inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
						wall_coords[wid].append(inpos)
				else:
					ncell=list(mesh.regions(2,wid))[0]
					centre=centroid(mesh,position,3,ncell)
					centre=[centre[0],centre[1]]
					w_points=[]
					for lid in mesh.borders(2,wid):
						for pid in mesh.borders(1,lid):
							w_points.append((position[pid][0],position[pid][1]))
					print 'w_points',w_points
					w_points=set(w_points)
					w_points=list(w_points)
					print 'w_points',w_points

					wall_coords[wid].append([w_points[0][0],w_points[0][1]])
			
					n=sqrt((w_points[0][0]-centre[0])*(w_points[0][0]-centre[0])+(w_points[0][1]-centre[1])*(w_points[0][1]-centre[1]))
					inpos=[centre[0]+(w_points[0][0]-centre[0])*(1-off/n),centre[1]+(w_points[0][1]-centre[1])*(1-off/n)]
					wall_coords[wid].append(inpos)
					n=sqrt((w_points[1][0]-centre[0])*(w_points[1][0]-centre[0])+(w_points[1][1]-centre[1])*(w_points[1][1]-centre[1]))
					inpos=[centre[0]+(w_points[1][0]-centre[0])*(1-off/n),centre[1]+(w_points[1][1]-centre[1])*(1-off/n)]
					wall_coords[wid].append(inpos)

					wall_coords[wid].append([w_points[1][0],w_points[1][1]])
	elif mesh.degree()==2:
		for pid,pos in db_position.iteritems():
			position[pid]=((pos[0]-min(xvals))/norm,(pos[1]-min(yvals))/norm)
		cell_to_coordinates={}
		for cid in mesh.wisps(2):
			cell_to_coordinates[cid]=[]
			for pid in ordered_pids(mesh,cid):
				cell_to_coordinates[cid].append([position[int(pid)][0],position[int(pid)][1]])
			

		wall=db.get_property('wall')
		edge_coords={}
		wall_coords={}
		off=0.01
		from math import atan2
		for eid in graph.edges():
			edge_coords[eid]=[]
			wid=wall[eid]
			centre=centroid(mesh,position,2,graph.source(eid))
			for pid in mesh.borders(1,wid):
				pidpos=[position[int(pid)][0],position[int(pid)][1]]

				edge_coords[eid].append(pidpos)

				ngbrs=list(mesh.region_neighbors(0,pid))
				pt0=position[ngbrs[0]]
				pt1=position[ngbrs[1]]
				angle0=atan2((pt0[1]-pidpos[1]),(pt0[0]-pidpos[0]))
				angle1=atan2((pt1[1]-pidpos[1]),(pt1[0]-pidpos[0]))
				adiff=angle0-angle1				
				
				#*
				n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
				inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
				#*
				edge_coords[eid].append(inpos)
			edge_coords[eid][3]=edge_coords[eid][2]
			edge_coords[eid][2]=inpos

			if wid in list(wall_coords.iterkeys()):
				wall_coords[wid].append(edge_coords[eid][3])
				wall_coords[wid].append(edge_coords[eid][2])
				wall_coords[wid].append(edge_coords[eid][1])
			else:
				wall_coords[wid]=[edge_coords[eid][0],edge_coords[eid][1],edge_coords[eid][2]]
		for wid in mesh.wisps(1):
			if wid not in list(wall_coords.iterkeys()):
				wall_coords[wid]=[]
				centre=centroid(mesh,position,2,list(mesh.regions(mesh.degree()-1,wid))[0])
				for pid in mesh.borders(1,wid):
					pidpos=[position[int(pid)][0],position[int(pid)][1]]
					wall_coords[wid].append(pidpos)
				
					#*
					n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
					inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
					#*
					wall_coords[wid].append(inpos)
				wall_coords[wid][3]=wall_coords[wid][2]
				wall_coords[wid][2]=inpos


	for v in range (0,max(list(vindex.itervalues()))+1):

		if len(prop_list)==11:
			fig=figure(figsize=(60,45)) 
			subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
		        wspace=0.2, hspace=0.2)
			spno=341
			text_coords={0:[0.05,0.926],1:[0.28,0.926],2:[0.51,0.926],3:[0.74,0.926],4:[0.05,0.6125],5:[0.28,0.6125],6:[0.51,0.6125],7:[0.74,0.6125],8:[0.05,0.29],9:[0.28,0.29],10:[0.51,0.29]}			
		if len(prop_list)==9:
			fig=figure(figsize=(60,45)) 
			subplots_adjust(left=0.13, bottom=0.05, right=0.87, top=0.95,
		        wspace=0.2, hspace=0.2)
			spno=331
			text_coords={0:[0.125,0.926],1:[0.39,0.926],2:[0.65,0.926],3:[0.125,0.6125],4:[0.39,0.6125],5:[0.65,0.6125],6:[0.125,0.29],7:[0.39,0.29],8:[0.65,0.29]}
			lwidth=8
			textsize=80
		if len(prop_list)==6:
			fig=figure(figsize=(60,45)) 
			subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.85,
		        wspace=0.2, hspace=0.2)
			spno=231
			text_coords={0:[0.05,0.825],1:[0.35,0.825],2:[0.675,0.825],3:[0.05,0.45],4:[0.35,0.45],5:[0.675,0.45]}
			lwidth=8
			textsize=80
		if len(prop_list)==5:
			fig=figure(figsize=(60,45)) 
			subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.85,
		        wspace=0.2, hspace=0.2)
			spno=231
			text_coords={0:[0.05,0.825],1:[0.35,0.825],2:[0.675,0.825],3:[0.05,0.45],4:[0.35,0.45]}
			lwidth=8
			textsize=80	
		if len(prop_list)==4:
			fig=figure(figsize=(60,45)) 
			subplots_adjust(left=0.1125, bottom=0.05, right=1-0.1125, top=0.95,
		        wspace=0.2, hspace=0.2)
			spno=221
			text_coords={0:[0.1,0.9],1:[0.5,0.9],2:[0.1,0.45],3:[0.5,0.45]}	
			lwidth=8
			textsize=80		
		if len(prop_list)==3:
			fig=figure(figsize=(45,13)) 
			spno=131
			text_coords={0:[0.1,0.9],1:[0.375,0.9],2:[0.65,0.9]}	
		if len(prop_list)==2:
			fig=figure(figsize=(60,45))
			subplots_adjust(left=0.05, bottom=0.25, right=0.95, top=0.75,
		        wspace=0.2, hspace=0.2)
			spno=121
			text_coords={0:[0.05,0.70],1:[0.55,0.70]}
			lwidth=8
			textsize=80
		if len(prop_list)==1:
			fig=figure(figsize=(30,30)) 
			spno=111
			text_coords={0:[0.1,0.9]}
			lwidth=8
			textsize=80

		pat = []
		for ncid,cid in enumerate(mesh.wisps(mesh.degree())):
			if vindex[cid]==v:
			    polygon = Polygon(cell_to_coordinates[cid], True)
			    polygon.set_linewidth(1.0)
			    pat.append(polygon)
		epat = []
		for ncid,eid in enumerate(graph.edges()):
			if vindex[graph.source(eid)]==v and vindex[graph.target(eid)]>=vindex[graph.source(eid)]:
				polygon = Polygon(edge_coords[eid], True)
				polygon.set_linewidth(1.0)
				epat.append(polygon)
		
		wpat = []
		for ncid,wid in enumerate(mesh.wisps(mesh.degree()-1)):
			app_wid=True
			r_cells = list(mesh.regions(mesh.degree()-1,wid))
			for cid in r_cells:
				if vindex[cid]!=v:
					app_wid=False
			if app_wid:
				polygon = Polygon(wall_coords[wid], True)
				polygon.set_linewidth(1.0)
				wpat.append(polygon)
		


		i=0
		species_desc=db.get_property('species_desc')
		print 'sd',species_desc
		for prop in prop_list:
			# auxin graph
			if i<9:
				ax=fig.add_subplot(spno)
			else:
				if i==9:
					ax=fig.add_subplot(3,4,10)
				if i==10:
					ax=fig.add_subplot(3,4,11)
			ax.set_xticks([])
			ax.set_yticks([])
			box(on=None)

			xlim(0,1)    # set the xlim to xmin, xmax
			ylim(0,1)    # set the ylim to ymin, ymax
			cprop=db.get_property(prop)
			Ax_colors = []
			if species_desc[prop]=='cell' or species_desc[prop]=='CELL':
				for ncid,cid in enumerate(mesh.wisps(mesh.degree())):
					if vindex[cid]==v:
				    		Ax_colors.append(cprop[cid])
				p = PatchCollection(pat, cmap=matplotlib.cm.jet, alpha=1.,linewidths=(lwidth,))
			elif species_desc[prop]=='edge' or species_desc[prop]=='EDGE':
				for eid in graph.edges():
					sid=graph.source(eid)
					tid=graph.target(eid)
					if vindex[sid]==v and vindex[sid]<=vindex[tid]:
						Ax_colors.append(cprop[eid])
				p = PatchCollection(epat, cmap=matplotlib.cm.jet, alpha=1.,linewidths=(5,))
				Lx_colors=[]
				for ncid,cid in enumerate(mesh.wisps(mesh.degree())):
					if vindex[cid]==v:
				    		Lx_colors.append(1)
				p0 = PatchCollection(pat, cmap=matplotlib.cm.gray, alpha=1.,linewidths=(5,))
				p0.set_array(array(Lx_colors))
				p0.set_clim(0,1)
				ax.add_collection(p0)
			
			elif species_desc[prop]=='wall' or species_desc[prop]=='WALL':
				for ncid,wid in enumerate(mesh.wisps(mesh.degree()-1)):
					app_wid=True
					r_cells = list(mesh.regions(mesh.degree()-1,wid))
					for cid in r_cells:
						if vindex[cid]!=v:
							app_wid=False
					if app_wid:
						Ax_colors.append(cprop[wid])
				p = PatchCollection(wpat, cmap=matplotlib.cm.jet, alpha=1.,linewidths=(5,))
			

			else:
				print prop,'type not CELL,EDGE or WALL'

			
			p.set_array(array(Ax_colors))
			p.set_clim(cmaps[i][0],cmaps[i][1])
			ax.add_collection(p)
			cb=colorbar(p,shrink=0.7)
			for t in cb.ax.get_yticklabels():
	     			t.set_fontsize(40)
			cb.set_ticks([cmaps[i][0],(cmaps[i][1]-cmaps[i][0])/4,(cmaps[i][1]-cmaps[i][0])/2,3*(cmaps[i][1]-cmaps[i][0])/4,cmaps[i][1]])
			cb.set_ticklabels([cmaps[i][0],(cmaps[i][1]-cmaps[i][0])/4,(cmaps[i][1]-cmaps[i][0])/2,3*(cmaps[i][1]-cmaps[i][0])/4,cmaps[i][1]])
			figtext(text_coords[i][0], text_coords[i][1], prop,size=textsize)
			spno += 1
			i += 1

		if mesh.degree()==3:
			if iteration >0:
				fname = '%s/cmap_%05d_%02d.png'%(savename,iteration,v)
			else:
				fname = '%s.png' % savename
		elif mesh.degree()==2:
			if iteration >0:
				fname = '%s/cmap_%05d.png'%(savename,iteration)
			else:
				fname = '%s.png' % savename
		savefig(fname)
		clf()
		close(fig)


