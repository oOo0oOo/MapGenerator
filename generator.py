from geometry import *
import matplotlib.pyplot as plt
from matplotlib import patches
from os import listdir
import pygame
import sys
import time
import random

class Map(object):
	def __init__(self, size_x, size_y):
		self.size_x = size_x
		self.size_y = size_y

		# Setup the playing field
		self.board = Rectangle(Point(size_x,size_y), Point(size_x * 2, size_y * 2))

	def generate_map(self, num_sectors, num_houses = 0, num_hedges = 0):
		# Create player HQs
		self.players = [self.board.random_point(0.7)]
		self.players.append(self.board.mirror_center(self.players[0]))

		#Add houses
		self.houses = []
		num = int(round(float(num_houses)/2))

		for i in range(num):
			h = self.board.random_point(max_dist_center = 0.85)
			self.houses.append(h)
			self.houses.append(self.board.mirror_center(h))

		#Add hedges
		self.hedges = []
		if num_hedges:
			max_dist = min([self.size_x, self.size_y])/num_hedges

			for i in range(int(round(float(num_hedges)/2))):
				num_points = random.randrange(2, 3)
				points = [self.board.random_point(max_dist_center = 0.85)]
				while len(points) < num_points:
					p = self.board.random_point(max_dist_center = 0.85)
					if Vector(p, points[len(points)-1]).length() < max_dist:
						points.append(p)

				self.hedges.append(points)
				self.hedges.append([self.board.mirror_center(p) for p in points])

		self.sectors = self.board.create_mesh(num_sectors, self.players)

	def draw_map_pygame(self):
		pygame.init() 
		window = pygame.display.set_mode((self.size_x, self.size_y))
		pygame.display.set_caption('Random Map for MapWarfare')

		transf_vector = Point(-self.size_x, -self.size_y)
		sector_font = pygame.font.SysFont("Tahoma", 35)

		def draw_bg():
			window.fill((15, 45, 20))

		def draw_image(img_path, point, rotation = 0, scale = 1):
			'''Rotation is around point supplied'''

			img = pygame.image.load(img_path)
			p = point.transform(transf_vector)

			rotated = pygame.transform.rotozoom(img, rotation, scale)
			#get the rect of the rotated surf and set it's center to the oldCenter
			rotRect = rotated.get_rect()
			d1, d2 = rotRect.size

			p.transform(Point(d1/2.0, d2/2.0))
			rotRect.center = (p.x, p.y)

			window.blit(rotated, rotRect)
			

		def draw_text(text, point):
			p = point.transform(transf_vector)
			label = sector_font.render(text, 50, (20,20,20))

			# Cetner on point
			rect = label.get_rect()
			d1, d2 = rect.size

			p.transform(Point(d1/2.0, d2/2.0))
			rect.center = (p.x, p.y)

			window.blit(label, rect)

		def draw_line(point1, point2, color = (0, 0, 0), width = 1):
			p1 = point1.transform(transf_vector)
			p2 = point2.transform(transf_vector)
			pygame.draw.line(window, color, (p1.x, p1.y), (p2.x, p2.y), width)

		draw_bg()

		# Draw borders
		for s in self.sectors:
			# Draw borders
			for l in s.sides:
				draw_line(l.p1, l.p2, (98, 98, 135), width = 3)

		# Hedges
		for points in self.hedges:
			last = len(points) - 1
			for i, p in enumerate(points):
				if i != last:
					draw_line(p, points[i+1], (80, 30, 0), 10)

		#Draw houses
		base = 'img/houses/'
		all_houses = [base + f for f in listdir(base) if f != 'hq.png']
		for p in self.houses:
			rot = random.randrange(0, 360)
			scale = random.randrange(65, 100)/100.0
			draw_image(random.choice(all_houses), p, rot, scale)


		# Draw players
		for p in self.players:
			draw_image('img/houses/hq.png', p, random.randrange(0, 360), 1.2)

		# Draw sector numbers
		for num, s in enumerate(self.sectors):
			# Draw number
			draw_text(str(num+1), s.centroid)


		'''
		#691F01
		'''

		pygame.display.flip()

		# input handling (somewhat boilerplate code):
		while True: 
		   for event in pygame.event.get(): 
		      if event.type == pygame.QUIT: 
		          sys.exit(0) 
		      else: 
		          pass

	def draw_map_matplotlib(self):
		#Create figure
		fig = plt.figure(1)

		def draw_point(point, color = 'b', marker = 'o', size = 20):
			plt.scatter([point.x], [point.y], color = color, marker = marker, s = size)

		def draw_vector(vector, color = 'k'):
			draw_line(vector.p1, vector.p2, color)

		def draw_line(p1, p2, color = 'k', linestyle = '-', linewidth = 2):
			plt.plot([p1.x, p2.x], [p1.y, p2.y], color=color, linestyle=linestyle, linewidth=linewidth)

		def draw_circle(center, radius, color = 'r', fig = fig):
			circle = plt.Circle((center.x, center.y), radius, color=color, fill=False)
			fig.gca().add_artist(circle)

		def draw_font(point, text):
			plt.text(point.x, point.y, text, size = 20, color = 'r')

		def draw_polygon(polygon, color = 'k', linestyle = '-', linewidth = 2):
			polygon.reorder_points()
			[draw_line(s.p1, s.p2, color = color, linestyle = linestyle, linewidth = linewidth) for s in polygon.sides]

		def draw_polygon_fill(polygon):
			polygon.reorder_points()
			pts = [p.coordinates() for p in polygon.points]
			poly = patches.Polygon(pts,facecolor=random.choice(['r', 'k', 'b', 'y', 'm', 'g', 'c']))
			plt.gca().add_patch(poly)

		def draw_rectangle(rectangle, color = 'k', linestyle = '-', linewidth = 2):
			[draw_line(s.p1, s.p2, color = color, linestyle = linestyle, linewidth = linewidth) for s in rectangle.sides]

		# Plot field
		draw_rectangle(self.board)

		# Plot sectors
		[draw_polygon(s, linestyle = '--', linewidth = 2) for s in self.sectors]

		# Plot sector number
		[draw_font(s.centroid, str(i+1)) for i, s in enumerate(self.sectors)]

		# Plot houses
		for h in self.houses:
			draw_point(h, color = 'k', marker = 's', size = random.randrange(120, 350))

		# Plot hedges
		for points in self.hedges:
			last = len(points) - 1
			for i, p in enumerate(points):
				if i != last:
					draw_line(p, points[i+1], '#691F01', linewidth = 9)

		# Plot players
		[draw_point(p, color = 'b', marker = 'D', size = 250) for p in self.players]

		normal = [self.size_x, self.size_x * 2, self.size_y, self.size_y * 2]
		large = [0, self.size_x * 3, 0, self.size_y * 3]
		plt.axis(normal, 'equal')

		plt.show()

if __name__ == '__main__':
	m = Map(1200, 600)
	m.generate_map(10, 10, 4)
	#m.draw_map_matplotlib()
	m.draw_map_pygame()
