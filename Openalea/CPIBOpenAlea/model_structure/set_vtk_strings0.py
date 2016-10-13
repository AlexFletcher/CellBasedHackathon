def set_vtk_strings0(db):

	cell_wall_width=db.get_property('cell_wall_width')

	mesh=get_mesh(db)

	if mesh.degree()==2:
		n_cells=len(list(mesh.wisps(mesh.degree())))

		type_string=""
		for i in range(n_cells):
			type_string +="7 "

		pid_trans={}
		position_string=""
		interior_points = interior_offset(db,cell_wall_width)
		str_pos=0
		for cid,pts in interior_points.iteritems():
			pid_trans[cid]=[]
			for pt in pts:
				position_string += "%s %s %s " % (str(pt[0]),str(pt[1]),str(pt[2]))
				pid_trans[cid].append(str_pos)
				str_pos+=1
		n_points=str_pos

		connect_string=""
		offset_string=""
		current_offset=0
		cid_trans={}
		cid_pos=0
		for cid,pts in interior_points.iteritems():
			current_offset +=len(pts)
			pidpos=0
			for pt in pts:
				connect_string += "%s " % str(pid_trans[cid][pidpos])
				pidpos+=1
			offset_string+="%s " % str(current_offset)
			cid_trans[cid_pos]=cid
			cid_pos+=1
	elif mesh.degree()==3:

		cells_walls = {}
		n_cells=0
		for cid in mesh.wisps(3):
			cells_walls[cid]=list(mesh.borders(3,cid))
			n_cells+=len(list((mesh.borders(3,cid))))


		template_db=db.get_property('template_db')
		interior_points = interior_offset(template_db,cell_wall_width)
		position = db.get_property('position')
		zpts={}
		for pid,pos in position.iteritems():
			if pos[2] not in zpts.iterkeys():
				zpts[pos[2]]=[pid]
			else:
				zpts[pos[2]].append(pid)


		pid_trans={}
		position_string=""
		str_pos=0
		cid_off={}
		zvals=sorted(list(zpts.iterkeys()))
		for cid,pts in interior_points.iteritems():
			pid_trans[cid]=[]
			cid_off[cid]=len(pts)
			for zpos in range(len(zvals)-1):
				for pt in pts:
					position_string += "%s %s %s " % (str(pt[0]),str(pt[1]),str(zvals[zpos]+cell_wall_width))
					pid_trans[cid].append(str_pos)
					str_pos+=1
				for pt in pts:
					position_string += "%s %s %s " % (str(pt[0]),str(pt[1]),str(zvals[zpos+1]-cell_wall_width))
					pid_trans[cid].append(str_pos)
					str_pos+=1
		n_points=str_pos


		connect_string=""
		offset_string=""
		current_offset=0
		cid_trans={}
		cid_pos=0
		offset_count=0
		for cid,pts in interior_points.iteritems():

			pidpos=0
			for zpos in range(len(zvals)-1):
				current_offset +=len(pts)					
				for pt in pts:
					connect_string += "%s " % str(pid_trans[cid][pidpos])
					pidpos+=1
				offset_string+="%s " % str(current_offset)
				offset_count +=1
				current_offset +=len(pts)
				for pt in pts:
					connect_string += "%s " % str(pid_trans[cid][pidpos])
					pidpos+=1
				offset_string+="%s " % str(current_offset)
				offset_count +=1
			cid_trans[cid_pos]=cid
			cid_pos+=1
		
		cpos=0
		for cid,pts in interior_points.iteritems():
			
			N=len(pts)
			Z=len((list(zpts.iterkeys())))-1
			for zpos in range(Z):
				for pt in range(N):
					connect_string += "%s " % str(cpos + pt%N)
					connect_string += "%s " % str(cpos+(pt+1)%N)
					connect_string += "%s " % str(cpos+N+(pt+1)%N)
					connect_string += "%s " % str(cpos+N+pt%N)
					current_offset += 4
					offset_string+="%s " % str(current_offset)
					offset_count +=1
				cpos+=2*N
			#cpos+=N
						
		n_cells=offset_count

		type_string=""
		for i in range(n_cells):
			type_string +="7 "


	db.set_property('vtk_strings',{'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,'n_cells':str(n_cells),'n_points':str(n_points)})
	db.set_description('vtk_strings','description of tissue for vtk format')
