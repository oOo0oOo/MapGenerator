import math
import random
import numpy as np
from scipy.spatial import Voronoi, ConvexHull

class Point(object):
	def __init__(self, x, y):
		self.x = float(x)
		self.y = float(y)

	def __str__(self):
		return 'Point: ({}, {})'.format(self.x, self.y)

	def coordinates(self):
		'''Returns list: [x, y]'''
		return [float(self.x), float(self.y)]

	def transform(self, vector, ratio = 1):
		x = self.x + ratio * vector.x
		y = self.y + ratio * vector.y
		return Point(x, y)

	def copy(self):
		return Point(float(self.x), float(self.y))

class Vector(object):
	def __init__(self, p1, p2):
		self.p1 = p1.copy()
		self.p2 = p2.copy()
		self.vect = Point(self.p2.x - self.p1.x, self.p2.y - self.p1.y)

	def get(self):
		return self.vect.copy()

	def length(self):
		return math.sqrt((self.vect.x**2) + (self.vect.y**2))

	def relative_point(self, ratio):
		if 0 <= ratio <= 1:
			x = float(self.p1.x + ratio * self.vect.x)
			y = float(self.p1.y + ratio * self.vect.y)
			return Point(x, y)

	def intersect(self, vector):
		'''
			Intersect another vector with this vector
		'''
		p3 = vector.p1
		p4 = vector.p2
		p1 = self.p1
		p2 = self.p2

		x = (p1.x*p2.y - p1.y*p2.x)*(p3.x - p4.x) - (p1.x - p2.x)*(p3.x*p4.y - p3.y*p4.x)
		x /= ((p1.x - p2.x)*(p3.y - p4.y) - (p1.y-p2.y)*(p3.x-p4.x))

		y = ((p1.x*p2.y - p1.y*p2.x)*(p3.y - p4.y)) - ((p1.y - p2.y)*(p3.x*p4.y - p3.y*p4.x))
		y /= ((p1.x - p2.x)*(p3.y - p4.y) - (p1.y-p2.y)*(p3.x-p4.x))

		return Point(x, y)

	def on_vector(self, point):
		'''
			Check if point is on vector. 
			point = p1 + u * vect
		'''
		uxx, uyy = True, True
		ux, uy = 0, 0
		if self.vect.x != 0:
			ux = (point.x - self.p1.x) / self.vect.x
		else:
			uxx = False

		if self.vect.y != 0:
			uy = (point.y - self.p1.y) / self.vect.y
		else:
			uyy = False

		if ux == uy and 0 <= ux <= 1 and uxx and uyy:
			return True
		elif uxx and 0 <= ux <= 1:
			return True
		elif uyy and 0 <= uy <= 1:
			return True
		else:
			return False

	def circle_collision(self, point, radius):
		minus = Point(self.p2.x - self.p1.x, self.p2.y - self.p1.y)
		diff = Point(self.p1.x - point.x, self.p1.y - point.y)

		a = minus.x **2 + minus.y **2
		b = 2 * ((minus.x * diff.x) + (minus.y * diff.y))
		c = (diff.x **2) + (diff.y **2) - (radius **2)
		delta = b **2 - (4 * a * c)
		
		points = []
		if delta < 0: # No intersection
			pass
		else:
			root = math.sqrt(delta)
			for i in [-1, 1]:
				t = (i * -b + root) / (2 * a)
				x = self.p1.x + (i * minus.x * t)
				y = self.p1.y + (i * minus.y * t)
				points.append(Point(x, y))

		return points

	def closest_point(self, point):
		u = ((point.x - self.p1.x) * (self.p2.x - self.p1.x)) + ((point.y - self.p1.y) * (self.p2.y - self.p1.y))
		u /= ((self.p2.x - self.p1.x) ** 2 + (self.p2.y - self.p1.y) ** 2)
		x = self.p1.x + (u * (self.p2.x - self.p1.x))
		y = self.p1.y + (u * (self.p2.y - self.p1.y))
		p = Point(x, y)
		if self.on_vector(p):
			return p
		else:
			l1 = Vector(self.p1, p).length()
			l2 = Vector(self.p2, p).length()
			if l1 < l2:
				return self.p1.copy()
			else:
				return self.p2.copy()
		
class Polygon(object):
	def __init__(self, points):
		self.points = [p.copy() for p in points]
		self.update()

	def copy(self):
		return Polygon([p.copy() for p in self.points])

	def add_point(self, point):
		self.points.append(point.copy())
		self.reorder_points()

	def intersect_polygon(self, polygon):
		'''
			Checks whether polygon intersects with other polygon.
		'''
		for s in polygon.sides:
			for s_own in self.sides:
				
				try:
					p = s.intersect(s_own)
					if s.on_vector(p) and s_own.on_vector(p):
						return True

				except ZeroDivisionError:
					return False

		return False

	def update(self):
		# Calculate centroid & area
		num_points = len(self.points)
		self.centroid = Point(0, 0)
		self.area = 0.0
		self.sides = []

		# From each point to the following point
		for i, p1 in enumerate(self.points):
			p2 = self.points[(i+1)%num_points]
			self.sides.append(Vector(p1, p2))

			# Calculate area
			a = p1.x*p2.y - p2.x*p1.y
			self.area += a

			# Transform centroid
			vect = Point((p1.x + p2.x)*a, 
						(p1.y + p2.y)*a)
			self.centroid = self.centroid.transform(vect)

		# Scaling
		if self.area:
			self.area *= 0.5
			self.centroid.x /= (6 * self.area)
			self.centroid.y /= (6 * self.area)
		else:
			self.centroid = self.sides[0].relative_point(0.5)

	def reorder_points(self):
		'''Sorts polygon for display'''
		pp = [p.coordinates() for p in self.points]
		cent = self.centroid.coordinates()
		pp.sort(key=lambda p: math.atan2(p[1]-cent[1],p[0]-cent[0]))
		self.points = [Point(p[0], p[1]) for p in pp]
		self.update()

	def borders_polygon(self, polygon):
		'''
			Checks whether the polygon borders another polygon.
			Returns a list containing the common sides.
		'''

	def merge(self, polygon):
		'''
			Merges the polygon with another polygon. The two new polygons must have
			at least one side in common to be merged.
		'''
		pol_points = [tuple(p.coordinates()) for p in polygon.points]
		points = [tuple(p.coordinates()) for p in self.points]

		num_common = 0
		for p in points:
			if p in pol_points:
				num_common += 1

		if num_common < 2:
			raise RuntimeError('Polygons do not have at least one common side.')

		all_points = set(pol_points + points)
		if len(all_points) >= 3:
			arr = np.array([list(l) for l in all_points])
			hull = ConvexHull(arr).points
			self.points = [Point(p[0], p[1]) for p in hull]
			self.reorder_points()
		else:
			pass
			#raise RuntimeError('Couldnt merge too little points')

class Rectangle(object):
	def __init__(self, bl, tr):
		self.bl = bl
		self.tr = tr
		self.diagonal = Vector(self.bl, self.tr)
		self.center = self.diagonal.relative_point(0.5)

		# Add all sides as vectors
		self.sides = []
		br = Point(self.tr.x, self.bl.y)
		tl = Point(self.bl.x, self.tr.y)
		self.sides.append(Vector(self.bl, br))
		self.sides.append(Vector(self.bl, tl))
		self.sides.append(Vector(self.tr, br))
		self.sides.append(Vector(self.tr, tl))

		self.area = self.sides[0].length() * self.sides[1].length()

		# Corners
		self.corners = [Point(p.x, p.y) for p in [self.bl, br, self.tr, tl]]

	def inside(self, point):
		if self.bl.x <= point.x <= self.tr.x:
			if self.bl.y <= point.y <= self.tr.y:
				return True
		return False

	def closest_corner(self, point):
		lengths = [(Vector(n, point).length(), n) for n in self.corners]
		return sorted(lengths)[0][1]

	def intersect_wall(self, point, radius):
		all_points = []
		for s in self.sides:
			points = s.circle_collision(point, radius)
			for p in points:
				if self.inside(p):
					all_points.append(p)

		return all_points

	def intersect_walls(self, vector):
		ps = []
		for s in self.sides:
			p = s.intersect(vector)
			if s.on_vector(p) and vector.on_vector(p):
				ps.append(p)
		return ps

	def dist_center(self, point, relative = False):
		dist = Vector(point, self.center).length()
		if not relative:
			return dist
		else:
			return dist / (self.diagonal.length()/2)

	def mirror_center(self, point):
		vect = Vector(point, self.center).get()
		return point.transform(vect, 2)

	def random_point(self, min_dist_center = 0, max_dist_center = 1, relative = True): 
		'''
			Relative decides whether *_dist_center is 
			relative (0 -> 1) or 
			absolute (0 -> length_diagonal)
		'''

		def rp():
			x = random.uniform(self.bl.x, self.tr.x)
			y = random.uniform(self.bl.y, self.tr.y)
			return Point(x, y)

		if not relative:
			len_ok = min_dist_center < (self.diagonal.length()/2) * 0.99
			len_ok = len_ok and (max_dist_center < (self.diagonal.length()/2) * 0.99)
		else:
			len_ok = min_dist_center / (self.diagonal.length()/2) < 0.99
			len_ok = len_ok and (max_dist_center / (self.diagonal.length()/2) < 0.99)

		if (min_dist_center >= 0 or max_dist_center !=1) and len_ok:
			while True:
				p = rp()
				if max_dist_center > self.dist_center(p, relative) > min_dist_center:
					return p
		else:
			return rp()

	def create_mesh(self, num_sectors, players = []):
		if num_sectors < 4:
			print 'num_sectors > 4!'
			return []

		# Calculate centers
		if num_sectors < 8:
			scale = 1
		elif num_sectors < 25:
			scale = 1.5
		else:
			scale = 2.5

		min_dist = (self.diagonal.length() / num_sectors) * scale

		# The evil double loop
		while True:
			tries = 0
			centers = [self.random_point(0, 0.9)]
			while len(centers) < num_sectors and not tries > 250:
				tries += 1
				n = self.random_point()
				lengths = [Vector(n, p).length() for p in centers]
				if min(lengths) > min_dist:
					centers.append(n)
					tries = 0

			if len(centers) == num_sectors:
				pl = []
				for p in players:
					lengths = [Vector(n, p).length() for n in centers]
					pl.append(sum(lengths)/num_sectors)

				if (min(pl) * 1.025) > max(pl):
					break

		# Calculate Voronoi 
		points = [p.coordinates() for p in centers]

		def frange(s, e, j):
			while s < e:
				yield s
				s += j

		# Add points around map (evenly spaced)
		s_x = self.tr.x - self.bl.x
		s_y = self.tr.y - self.bl.x

		dist = (s_x + s_y) / (4 * num_sectors)

		for i in frange(self.bl.x, self.tr.x, dist):
			points.append([i, self.bl.y*0.5])
			points.append([i, self.tr.y*1.5])

		for i in frange(self.bl.y, self.tr.y, dist):
			points.append([self.bl.x*0.5, i])
			points.append([self.tr.x*1.5, i])

		vor = Voronoi(np.array(points))

		# Create list of polygons (sectors)
		sectors = []
		vert = vor.vertices

		for points in vor.regions:
			if len(points) >= 3:
				col = []
				ps = [Point(vert[p][0], vert[p][1]) for p in points if p != -1]				
				poly = Polygon(ps)
				poly.reorder_points()

				# Check if sector is on the map and if it has to be cut
				ins = [self.inside(p) for p in ps]
				#side_inter = [self.intersect_walls(s) for s in poly.sides]
				#print ins, side_inter
				'''
				if any(side_inter) and any(ins):
					# Cut polygon at sides
					# All points already inside the map are accepted
					points = [p for i, p in enumerate(poly.points) if ins[i]]

					# Add all intersection points
					for res in side_inter:
						if res:
							points.append(res[0])

					#print points
					poly2 = Polygon(points)
					poly2.reorder_points()
					sectors.append(poly2)

				elif all(ins):
					sectors.append(poly)

				'''	
				if any(ins):
					sectors.append(poly)

				'''
				new_ps = []
				for side in poly.sides:

					ps = self.intersect_walls(side)
					if len(ps) == 1:
						if self.inside(side.p1):
							new_ps.append(side.p1)
						else:
							new_ps.append(side.p2)
						new_ps.append(ps[0])
					elif len(ps) == 2:
						[new_ps.append(p) for p in ps]
					elif not ps and self.inside(side.p1):
						new_ps.append(side.p1)
						new_ps.append(side.p2)
				'''

				#ps = []
				#for p in new_ps:
				#	if p not in ps:
				#		ps.append(p)
				'''
				if len(new_ps) >= 3:
					sectors.append(Polygon(new_ps))
				'''
				

		# Merge small sectors
		'''
		to_del = []
		for i, s in enumerate(sectors):
			if s.area > ((0.01/num_sectors)*self.area):
				pass
			else:
				to_del.append(i)

		for i in to_del:
			s = sectors[i]
			dists =	[tuple([Vector(s.centroid, c.centroid).length(), j]) for j, c in enumerate(sectors) if j not in to_del]
			arranged = sorted(dists)[:3]
			for _, closest in arranged:
				try:
					sectors[closest].merge(s)
					break
				except RuntimeError:
					pass

		for index in sorted(to_del, reverse=True):
			del sectors[index]

		'''
		# Add corners to closest sector
		'''
		for c in self.corners:
			dists =	[tuple([Vector(s.centroid, c).length(), i]) for i, s in enumerate(sectors)]
			closest = sorted(dists)
			i = 0
			found = False
			while not found and i < len(sectors):
				ind = closest[i][1]
				s = sectors[ind]
				# Make copy in case of intersecting
				orig_sector = s.copy()

				# Add point
				s.add_point(c)
				
				for sect in sectors:
					if s.intersect_polygon(sect):
						sectors[ind] = orig_sector
						break
				else:
					found = True
				i += 1
		
		'''
		return sectors