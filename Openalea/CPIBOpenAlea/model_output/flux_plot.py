import os, sys
import cPickle as pickle
from sa_oa.tissueshape import centroid
from sa_oa.container import ordered_pids
from model_utils.celltissue_util import *

from model_utils.db_utilities import get_mesh,get_graph,get_wall_decomp

from numpy import *
from pylab import *
from math import atan2

import matplotlib
from matplotlib.patches import Polygon,FancyArrowPatch
from matplotlib.collections import PatchCollection,PolyCollection
import matplotlib.pyplot as plt


from math import atan, cos,pi,sin,sqrt
from datetime import datetime
import time
a=datetime.now()

from sa_vp.plantgl.ext.color import JetMap, GrayMap, GreenMap

def flux_plot(db,iteration,dir_name,maxlength):
	mesh=get_mesh(db)
	graph=get_graph(db)
	db_position = db.get_property('position')
	poslist=list(db_position.itervalues())
	xvals,yvals = [[z[i] for z in poslist] for i in (0,1)]
	norm=max(max(xvals)-min(xvals),max(yvals)-min(yvals))
	position={}
	for pid,pos in db_position.iteritems():
		position[pid]=(100*(pos[0]-min(xvals))/norm,100*(pos[1]-min(yvals))/norm)
	vindex = db.get_property('vindex')
	wall=db.get_property('wall')
	if mesh.degree()==3:
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
		
		off=10
		for eid in graph.edges():
			edge_coords[eid]=[]
			
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
			else:
				centre=centroid(mesh,position,2,wall[eid])
				sid=graph.source(eid)
				for pidpos in cell_to_coordinates[sid]:
					n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
					inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
					edge_coords[eid].append(inpos)	
	elif mesh.degree()==2:
		cell_to_coordinates={}
		for cid in mesh.wisps(2):
			cell_to_coordinates[cid]=[]
			for pid in ordered_pids(mesh,cid):
				cell_to_coordinates[cid].append([position[int(pid)][0],position[int(pid)][1]])
			

		wall=db.get_property('wall')
		edge_coords={}

		off=10
		for eid in graph.edges():
			edge_coords[eid]=[]
			wid=wall[eid]
			centre=centroid(mesh,position,2,graph.source(eid))
			for pid in mesh.borders(1,wid):
				pidpos=[position[int(pid)][0],position[int(pid)][1]]
				edge_coords[eid].append(pidpos)
				n=sqrt((pidpos[0]-centre[0])*(pidpos[0]-centre[0])+(pidpos[1]-centre[1])*(pidpos[1]-centre[1]))
				inpos=[centre[0]+(pidpos[0]-centre[0])*(1-off/n),centre[1]+(pidpos[1]-centre[1])*(1-off/n)]
				edge_coords[eid].append(inpos)
			edge_coords[eid][3]=edge_coords[eid][2]
			edge_coords[eid][2]=inpos

	auxin_flux=db.get_property('PIN')
	fluxes=list(auxin_flux.itervalues())
	abs_fluxes=[]
	for i in fluxes:
		abs_fluxes.append(abs(i))
	max_flux=max(abs_fluxes)
	S=db.get_property('S')
	for eid,flux in auxin_flux.iteritems():
		if abs(flux)==max_flux:
			maxlength=S[wall[eid]]
	
	for v in range (0,max(list(vindex.itervalues()))+1):

		fig=figure(figsize=(30,30)) 
		spno=111
		text_coords={0:[0.1,0.9]}
	
	

		pat = []
		for ncid,cid in enumerate(mesh.wisps(mesh.degree())):
			if vindex[cid]==v:
			    polygon = Polygon(cell_to_coordinates[cid], True)
			    polygon.set_linewidth(100.0)
			    pat.append(polygon)
		arrows = []
		wall_decomp=get_wall_decomp(db)
		
		for wid,(eid1,eid2) in wall_decomp.iteritems():
			if vindex[graph.source(eid1)]==v and vindex[graph.target(eid2)]==v:
				pids=list(mesh.borders(1,wid))
				flux=auxin_flux[eid1]


				poss=[tuple(position[pids[0]]),tuple(position[pids[1]])]
				
				if poss[0][0]>poss[1][0]:
					x2,y2=poss[0][0],poss[0][1]
					x1,y1=poss[1][0],poss[1][1]
					angle=atan2(y2-y1,x2-x1)
				elif poss[0][0]<poss[1][0]:
					x1,y1=poss[0][0],poss[0][1]
					x2,y2=poss[1][0],poss[1][1]
					angle=atan2(y2-y1,x2-x1)
				elif poss[0][0]==poss[1][0]:
					x1,y1=poss[0][0],poss[0][1]
					x2,y2=poss[1][0],poss[1][1]
					angle=pi/2
				mp = (x1+(x2-x1)/2,y1+(y2-y1)/2)
				if max_flux!=0.0:
					alength=maxlength*sqrt(abs(flux)/max_flux)

					apoint1=(mp[0]-alength*sin(angle),mp[1]+alength*cos(angle))
					apoint2=(mp[0]+alength*sin(angle),mp[1]-alength*cos(angle))

					if flux>0:
						sc=centroid(mesh,position,2,graph.source(eid2))
					else:
						sc=centroid(mesh,position,2,graph.source(eid1))

					sdist1=sqrt((apoint1[0]-sc[0])**2+(apoint1[1]-sc[1])**2)
					sdist2=sqrt((apoint2[0]-sc[0])**2+(apoint2[1]-sc[1])**2)

					pt2=(mp[0]+alength*cos(angle)/2,mp[1]+alength*sin(angle)/2)
					pt3=(mp[0]-alength*cos(angle)/2,mp[1]-alength*sin(angle)/2)

					if sdist2>sdist1:
											
						arrows.append(Polygon([apoint2,pt2,pt3],True))
					else:			
						arrows.append(Polygon([apoint1,pt2,pt3],True))
		poslist=list(position.itervalues())
		xvals,yvals = [[z[i] for z in poslist] for i in (0,1)]


		
		species_desc=db.get_property('species_desc')
		ax=fig.add_subplot(spno)
		
		ax.set_xticks([])
		ax.set_yticks([])
		box(on=None)

		xlim(0,100 )    # set the xlim to xmin, xmax
		ylim(0, 100 )    # set the ylim to ymin, ymax
		#cprop=db.get_property(prop)
		Ax_colors = []
		for ncid,cid in enumerate(mesh.wisps(mesh.degree())):
			Ax_colors.append(1.0)
			p0 = PatchCollection(pat, cmap=matplotlib.cm.gray, alpha=1.,linewidths=(1,))
		
		
		
		
		p0.set_array(array(Ax_colors))
		p0.set_clim(0,1)
		ax.add_collection(p0)

		Flux_colours = []
		for wid,(eid1,eid2) in wall_decomp.iteritems():
			flux=auxin_flux[eid1]
			if vindex[graph.source(eid1)]==vindex[graph.source(eid2)]:
				if flux>0:
					Flux_colours.append(abs(auxin_flux[eid2]))
				else:
					Flux_colours.append(abs(auxin_flux[eid1]))

		
		if len(arrows)>0:
			p = PatchCollection(arrows, cmap=matplotlib.cm.jet, alpha=1.,linewidths=(1,))
			p.set_clim(0.0,max_flux)
			p.set_array(array(Flux_colours))
			ax.add_collection(p)
		
		#cb=colorbar(p0,shrink=0.7,anchor=(0.0,-0.5))


		if mesh.degree()==3:
			fname = '%s/flux_%05d_%02d.png'%(dir_name,iteration,v)
		elif mesh.degree()==2:
			fname = '%s/flux_%05d.png'%(dir_name,iteration)
		savefig(fname)
		clf()
		close(fig)


