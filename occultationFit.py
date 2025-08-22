# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 12:58:29 2024

@author: antoine
"""

from astroquery.vizier import Vizier
import xml.etree.ElementTree as ET
import glob
import os
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.coordinates import Longitude
import numpy as np
from astroquery.jplhorizons import Horizons
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import mpl_toolkits.mplot3d.art3d as art3d
from PIL import Image
from mpl_toolkits.mplot3d.art3d import PolyCollection
from scipy.optimize import minimize
import argparse
import sys
import re
import copy
import time
import random

np.seterr(all='ignore')

plt.ion()

file_path = os.path.dirname(os.path.realpath(__file__)) + '/'

class Event:
    
    def __init__(self, path, shouldRequestStarPosToDatabase = False):
        self.path = path 
        self.tree = ET.parse(path)
        self.root = self.tree.getroot()
        self.date = -1
        self.star = {}
        self.ast = {}
        self.observer = {}
        self.shouldRequestStarPosToDatabase = shouldRequestStarPosToDatabase
        
        self.__parseEvent()
    
        self.index = 0

    def getAst(self):
        return '(' + str(self.ast["id"]) + ')' + ' ' + self.ast["name"]  
    
    def __getitem__(self, observer_id):
        
        if observer_id not in self.observer:
            print("ID not found")
            return None
        
        return self.observer[observer_id]
    
    def __len__(self):
        return len(self.observer)
    
    def __iter__(self):
       return self

    def __next__(self):
        
       if self.index < len(self.observer):
                      
           keys = list(self.observer.keys())
     
           item = self[keys[self.index]]
           self.index += 1
         
           return item
       else:
      
           self.index = 0
     
           raise StopIteration
           
    
    def __parseEvent(self):

        self.__parseDate()
        self.__parseStar()
        self.__parseAst()
        self.__parseObserver()
    
    def __floatingTimeToTime(self, time):
        
        dec = time % 1        

        h = int(time - dec)
        m = int(dec * 60 - (dec * 60) % 1)
        s = ((dec * 60) % 1 ) * 60
        
        return h, m, s
    
    def __parseDate(self):
        
    
        d = self.findBalise(self.root, "Date")[0].text
        list_delement = d.split('|')
        d = '-'.join(list_delement[:-1])
        
        time = float(list_delement[-1])

        h, m, s = self.__floatingTimeToTime(time)
        
        d = str(d) + " " + str(h) + ":" + str(m) + ":" + str(s)
        
        self.date = Time(d, format = "iso")
     
    def __parseStar(self):
        main = self.findBalise(self.root, "Star")[0].text.split("|")
        
        self.star = {"catName" : main[0],
                     "catNum" : main[1],
                     "gaia_id" : main[3],
                     "ra" : float(main[4]) * 15,
                     "dec" : float(main[5])}
       
        
        if self.shouldRequestStarPosToDatabase:
            
            print("query star position to database")
            
            if self.star['catName'].lower() == "ucac4":
                ra, dec = self.__get_star_coordinates("ucac4", self.star["catNum"])
            elif self.star["catName"].lower() == "tycho2":
                ra, dec =self.__get_star_coordinates("tycho2", self.star["catNum"])
            elif self.star["catName"].lower() == "hip":
                ra, dec =self.__get_star_coordinates("hipparcos", self.star["catNum"])
            else:
                ra, dec =self.__get_star_coordinates("gaia", self.star["gaia_id"])
                
            self.star["ra"] = ra
            self.star["dec"] = dec
     
    def __parseAst(self):
        main = self.findBalise(self.root, "Asteroid")[0].text.split("|")
        self.ast = {"id" : int(main[0]),
                    "name" : main[1]}
        
    def __parseObserver(self):
        
        obss = self.findBalise(self.root, "Observer") 
        
        for obs in obss:
             
            list_of_element = self.findBalise(obs, "ID")[0].text.split("|")
            disa_list = self.findBalise(obs, "D")[0].text.split('|')
            reap_list = self.findBalise(obs, "R")[0].text.split("|")
            
            
            d_time_split = [x for x in disa_list[0].split(" ") if x != ""]
            r_time_split = [x for x in reap_list[0].split(" ") if x != ""]
            isPositif = [x for x in disa_list[1].split(' ') if x != ""]
            
            if "".join(isPositif) == "M":
                isPositif = False
                status = 1
            else:
                isPositif = True
                status = 0
            
            d_time = ":".join(d_time_split)
            d_err = 0
            
            r_time = ":".join(r_time_split)
            r_err = 0
            
            if len(r_time) == 1 or r_time == "." or d_time == ".":
                
                if r_time == ".":
                    r_time = None
                else:
                    d_time = None
                
                if isPositif:
                    status = 2    
            
                
            
            if not isinstance(r_time, type(None)) and r_time[-1] == '.':
                r_time += '0'
            if not isinstance(d_time, type(None)) and d_time[-1] == '.':
                d_time += '0'
            
            
                
            if len(disa_list[2]) != 0:
                d_err = float(disa_list[2])
            if len(reap_list[2]) != 0:
                r_err = float(reap_list[2])


            if not isinstance(d_time, type(None)):
                d_time = Time(self.date.value.split(" ")[0] + ' ' + d_time, format = "iso")
            else:
                d_time = None
                
            if not isinstance(r_time, type(None)):
                r_time = Time(self.date.value.split(" ")[0] + ' ' + r_time, format = "iso")
            else:
                r_time = None
                

            di = {"name1" : list_of_element[1], 
                  "name2" : list_of_element[2],
                  "long" : list_of_element[6],
                  "lat" : list_of_element[7],
                  "alt" : float(list_of_element[8]),
                  "isPositif" : isPositif, # doit être changer et supprimer au profit de status
                  "d" : d_time,
                  "d_err" : TimeDelta(d_err, format='sec'),
                  "r" : r_time,
                  "r_err": TimeDelta(r_err, format='sec'),
                  "status": status}
            
            self.observer[int(list_of_element[0])] = di
    
    def isDmissing(di):
        return di["status"] == 2 and isinstance(di["d"], type(None))
    
    def isRmissing(di):
        return di["status"] == 2 and isinstance(di["r"], type(None))
    
    def __get_star_coordinates(self, catalog_name, star_ids):
        # Initialiser un objet Vizier pour interroger les catalogues
        
        
        # Choisir le bon catalogue en fonction de son nom
        if catalog_name.lower() == 'tycho2':
            cat, const = self.__get_star_coordinates_tycho2(star_ids)
        elif catalog_name.lower() == 'hipparcos':
            cat, const = self.__get_star_coordinates_hip(star_ids)
        elif catalog_name.lower() == 'ucac4':
            cat, const = self.__get_star_coordinates_ucac(star_ids)
        else:
            cat, const = self.__get_star_coordinates_gaia(star_ids)
        
        tab = self.__vizier_request(cat, const)[0]
        
        return tab['_RAJ2000'][0], tab['_DEJ2000'][0]

    def __vizier_request(self, cat, const):
        v = Vizier(columns=['_RAJ2000', '_DEJ2000'], row_limit=10000)
        return v.query_constraints(catalog=cat, **const)


    def __get_star_coordinates_ucac(self, UCAC_id):
       
        catalog = 'I/322A'
        query_constraints = {
            'UCAC4': UCAC_id,
        } 

        return catalog, query_constraints

    def __get_star_coordinates_hip(self, hip_id):
       
        catalog = 'I/239'
        query_constraints = {
            'HIP': hip_id,
        } 

        return catalog, query_constraints

    def __get_star_coordinates_gaia(self, gaia_id):

        catalog = 'I/355'
        query_constraints = {
            'source_id': gaia_id,
        } 
        
        return catalog, query_constraints


    def __get_star_coordinates_tycho2(self, tycho_id):



        catalog = 'I/259'

        tyc_parts = tycho_id.split('-')
        tyc1 = tyc_parts[0]
        tyc2 = str(int(tyc_parts[1])) 
        tyc3 = tyc_parts[2]


        query_constraints = {
            'TYC1': tyc1,
            'TYC2': tyc2,
            'TYC3': tyc3
        } 

        return catalog, query_constraints
    
    
    def getId(self, ev):
        for i in self.observer:
            if ev == self.observer[i]:
                return i
            
    
    
    
    def findBalise(self, root, baliseTag):
        
        if root.tag == baliseTag and root != self.root:
            return root
        
        childs = []
        
        for child in root:
            c = self.findBalise(child, baliseTag)
            
            if not isinstance(c, list):
                c = [c]
            
            if len(c) != 0:
                childs = childs + c
                
        return childs
        
    
    def first(self):
        
        for key in self.observer:
            return self.observer[key]
    
    def findFirstPositif(self):
        
        for key in self.observer:
            if self.observer[key]["isPositif"]:
                return self.observer[key]

class ModelParser:
    def __init__(self, path, solution = None):
        
        self.solution = solution
        
        self.path = path
        
        self.nVertices = 0
        self.nTriangles = 0
        
        self.vertices = None
        self.vertices_transformed = None
        self.angle = 0
        
        self.triangles = []
        
        self.inertia_tensor = None
        self.center_of_mass = None
        self.volume = None
        self.principal_moments = None
        self.principal_axes = None
        
        self.id = random.randrange(1,1000)
        
        self.__parse()

    def __removeEmptyStrFromList(self, li):
        return [i for i in li if i != ""]
        
    def __parse(self):
       with open(self.path, 'r') as file:
           content = file.read()
       
       lines = content.split("\n")
      
       nVT = self.__removeEmptyStrFromList(lines[0].split(" "))
       
       self.nVertices = int(nVT[0])
       self.nTriangles = int(nVT[1])
       
       self.vertices = np.zeros((self.nVertices, 3))
       
       for i, line in enumerate(lines[1: self.nVertices + 1]):
           line = self.__removeEmptyStrFromList(line.split(" "))
           self.vertices[i, 0] = float(line[0])
           self.vertices[i, 1] = float(line[1])
           self.vertices[i, 2] = float(line[2])
       
       # Normalization step similar to MATLAB code
       norms = np.linalg.norm(self.vertices, axis=1)
       max_norm = max(norms) if max(norms) != 0 else 1.0
       self.vertices = self.vertices / max_norm
       
       # Parsing triangle indices, skipping every other line containing '3'
       i = self.nVertices + 1
       while i < len(lines):
           line = lines[i].strip()
           if line == "3":
               i += 1  # Skip this line
               continue
           elif line:
               cont = self.__removeEmptyStrFromList(line.split(' '))
               p1 = int(cont[0]) - 1  # Adjusting to 0-based index
               p2 = int(cont[1]) - 1
               p3 = int(cont[2]) - 1
               self.triangles.append([p1, p2, p3])
           i += 1
       
       self.triangles = np.array(self.triangles)
       self.compute_principal_moments()
    
    def angOfRot(self, time_of_interest):
        
        if self.solution == None:
            return None
        
        t_zero = self.solution["zero"]
        dt = (time_of_interest - t_zero) * 24 #jd
        n_period = dt / self.solution["period"]
        
        return 360 * (n_period - int(n_period)) 
        
  
    
    def rotate_along_spin(self, angle, onTransformed = False):
        
        angle = angle 
        
        da = self.angle - angle 
        self.angle = angle
   
        
        angle = np.deg2rad(da)
       
        mat = self.zRmat(angle)

        if onTransformed:
            self.vertices_transformed = self.vertices_transformed @ mat
        else:
            self.vertices = self.vertices @ mat
            
    def zRmat(self, angle):
        cos = np.cos(angle)
        sin = np.sin(angle)
        
        return np.asarray([[cos, -sin, 0], 
                          [sin,  cos, 0], 
                          [0,     0,  1]])
    
    def yRmat(self, angle):
        cos = np.cos(angle)
        sin = np.sin(angle)
        
        return np.asarray([[cos, 0, sin], 
                          [0,  1, 0], 
                          [-sin,     0,  cos]])
        

    def in_ecliptic_cord(self, lamb, beta, t, period, onTransformed = False):
        
        if onTransformed:
            vert = self.vertices_transformed.T
        else:
            vert = self.vertices.T

            
        lamb = np.deg2rad(lamb)        
        phi = np.deg2rad(float(self.solution["phase"]))
        t0 = float(self.solution["zero"])
        
        dt = (t - t0) * 24
        n_period = dt / self.solution["period"]
        

        beta =  np.deg2rad(beta)
        vert = self.zRmat(phi + np.deg2rad(360 * (n_period - int(n_period)))) @ vert
        vert = self.yRmat(np.pi/2 - beta) @ vert
        vert = self.zRmat(lamb) @ vert
        
        return vert.T
    
            
    def compute_volume_and_center_of_mass(self):
       volume = 0.0
       center_of_mass = np.zeros(3)
       
       for triangle in self.triangles:
           p1, p2, p3 = self.vertices[triangle]
           # Calculate vectors G and H for the triangle edges
           G = p2 - p1
           H = p3 - p1
           # Calculate normal N using cross product
           N = np.cross(G, H)
           # Triangle area
           dS = np.linalg.norm(N) / 2.0
           # Volume contribution (using scalar triple product)
           dV = np.dot(p1, N) / 6.0
           # Ensure volume contribution is positive for consistency with MATLAB
           dV = abs(dV)
           volume += dV
           # Calculate centroid of the triangle
           centroid = (p1 + p2 + p3) / 4.0  # Using (D + E + F) / 4 as in MATLAB
           # Update center of mass
           center_of_mass += dV * centroid
       
       self.volume = abs(volume)  # Ensure volume is positive
       if volume != 0:
           self.center_of_mass = center_of_mass / volume
       else:
           self.center_of_mass = center_of_mass
       
    def compute_inertia_tensor(self):
       if self.center_of_mass is None:
           self.compute_volume_and_center_of_mass()
       
       Ixx, Iyy, Izz, Ixy, Ixz, Iyz = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
       
       for triangle in self.triangles:
           p1, p2, p3 = self.vertices[triangle]
           # Translate vertices to the center of mass
           D = p1 - self.center_of_mass
           E = p2 - self.center_of_mass
           F = p3 - self.center_of_mass
           # Calculate vectors G and H for the triangle edges
           G = E - D
           H = F - D
           # Calculate normal N using cross product
           N = np.cross(G, H)
           # Volume contribution (using scalar triple product)
           dV = np.dot(D, N) / 6.0
           dV = abs(dV)  # Ensure volume contribution is positive for correct scaling
           
           # Moment contributions
           Dx, Dy, Dz = D
           Ex, Ey, Ez = E 
           Fx, Fy, Fz = F
           
           pom = dV / 20.0  # Use the same scaling factor as in MATLAB (dV / 20)
           dPxx = pom * (2 * Dx**2 + 2 * Ex**2 + 2 * Fx**2 + 2 * Dx * Ex + 2 * Dx * Fx + 2 * Ex * Fx)
           dPyy = pom * (2 * Dy**2 + 2 * Ey**2 + 2 * Fy**2 + 2 * Dy * Ey + 2 * Dy * Fy + 2 * Ey * Fy)
           dPzz = pom * (2 * Dz**2 + 2 * Ez**2 + 2 * Fz**2 + 2 * Dz * Ez + 2 * Dz * Fz + 2 * Ez * Fz)
           dPxy = pom * (2 * Dx * Dy + 2 * Ex * Ey + 2 * Fx * Fy + Dx * Ey + Dy * Ex + Dx * Fy + Dy * Fx + Ex * Fy + Ey * Fx)
           dPxz = pom * (2 * Dx * Dz + 2 * Ex * Ez + 2 * Fx * Fz + Dx * Ez + Dz * Ex + Dx * Fz + Dz * Fx + Ex * Fz + Ez * Fx)
           dPyz = pom * (2 * Dy * Dz + 2 * Ey * Ez + 2 * Fy * Fz + Dy * Ez + Dz * Ey + Dy * Fz + Dz * Fy + Ey * Fz + Ez * Fy)
           
           Ixx += dPyy + dPzz
           Iyy += dPxx + dPzz
           Izz += dPxx + dPyy
           Ixy -= dPxy
           Ixz -= dPxz
           Iyz -= dPyz
       
       # Inertia tensor
       self.inertia_tensor = np.array([[Ixx, Ixy, Ixz],
                                       [Ixy, Iyy, Iyz],
                                       [Ixz, Iyz, Izz]])
       
    def __organiseEigens(self, eigenVec, eigenVal):
        
        idxMax = []
        
        for i in range(eigenVec.shape[1]): # détermine ou sont les axes principaux et les inverse si ils sont négatif
            v = eigenVec[:, i]
            idxMax.append( np.where(np.abs(v) == np.max(np.abs(v)))[0][0] ) 
            
            if v[idxMax[-1]] == -1 * np.abs(v)[idxMax[-1]]:
                eigenVec[:, i] *= -1
            
        for i in range(len(idxMax)): #réordonne les vecteur ainsi que les valeur propre dans la diagonal
            

        
            if i != idxMax[i]:

                idxi = np.where(np.asarray(idxMax) == i)[0][0]
                tpmVec = eigenVec[:, idxi].copy()
                tpmVal = eigenVal[idxi].copy()
                tpmMax = idxMax[i]
                
                
                eigenVec[:, idxi] = eigenVec[:, i].copy()
                eigenVec[:, i] = tpmVec.copy()

                eigenVal[idxi] =  eigenVal[i].copy()
                eigenVal[i] = tpmVal.copy()
                
                idxMax[i] = i
                idxMax[idxi] = tpmMax
                
 
      
        return eigenVec, eigenVal                
   
    def compute_principal_moments(self, whatForm = 0): #whatForm représente la façon d'on les eigenvectors doivent être présenter afin d'imiter le comportement d'autre code déjà existant
                                                   #si 0 - la fonction vas faire en sorte que les axes principaux soit dans la diagonal et soit possitif
                                                   #si 1 - la fonction organisera les vecteurs propre pour qu'ils soit identique au script drawshape.m
       if self.inertia_tensor is None:
           self.compute_inertia_tensor()
       
       # Après calcul des vecteurs propres
       eigenvalues, eigenvectors = np.linalg.eig(self.inertia_tensor)

       if whatForm == 0:
           eigenvectors, eigenvalues = self.__organiseEigens(eigenvectors, eigenvalues) 
       elif whatForm == 1:
           idx = np.argsort(eigenvalues)
           eigenvalues = eigenvalues[idx]
           eigenvectors = eigenvectors[:, idx]
           eigenvectors[:, 2] *= -1
           eigenvectors[:, 1] *= -1
       
       # Stocker les moments et axes principaux
       self.principal_moments = eigenvalues
       self.principal_axes = eigenvectors
       
       # Transformer les sommets selon les axes principaux
       self.vertices_transformed = np.dot(eigenvectors, self.vertices.T).T
       self.angle = 0
      
    def plot3_vertices(self, transformed=False, asMatlab = True, block = False):
       fig = plt.figure(figsize=(15, 5))
       # Set the camera views similar to MATLAB code, with corrected order
       views = [(0, 0), (0, 90), (90, 0)]  # Corrected order: YZ, XY, XZ planes
       
       
       if asMatlab:
           self.compute_principal_moments(1)
       
       for i, view in enumerate(views):
           ax = fig.add_subplot(1, 3, i + 1, projection='3d')
           
           ax.set_facecolor('black')
           fig.patch.set_facecolor('black')
           
           vertices_to_plot = self.vertices_transformed if transformed else self.vertices
           tris = self.triangles
           
           # Create the plot with a grey surface and no lines
           ax.plot_trisurf(vertices_to_plot[:, 0], vertices_to_plot[:, 1], vertices_to_plot[:, 2], triangles=tris, 
                           color='grey', edgecolor='none', alpha=0.7)
           
           # num_triangles_to_exclude = min(10, len(tris))
           # tris_grey = tris[num_triangles_to_exclude:]  # Triangles à dessiner en gris (à partir du 11e)
           # tris_blue = tris[:num_triangles_to_exclude]  # Les 10 premiers triangles à dessiner en bleu
        
           # # Tracer le modèle gris sans les 10 premiers triangles
           # ax.plot_trisurf(vertices_to_plot[:, 0], vertices_to_plot[:, 1], vertices_to_plot[:, 2],
           #                 triangles=tris_grey, color='grey', edgecolor='none', alpha=0.8)
        
           # # Tracer les 10 premiers triangles en bleu
           # ax.plot_trisurf(vertices_to_plot[:, 0], vertices_to_plot[:, 1], vertices_to_plot[:, 2],
           #                 triangles=tris_blue, color='blue', edgecolor='none', alpha=1.0)
 
           
           # Add XYZ axes for better orientation
           max_range = np.array([vertices_to_plot[:, 0].max() - vertices_to_plot[:, 0].min(), 
                                 vertices_to_plot[:, 1].max() - vertices_to_plot[:, 1].min(), 
                                 vertices_to_plot[:, 2].max() - vertices_to_plot[:, 2].min()]).max() / 2.0
           mid_x = (vertices_to_plot[:, 0].max() + vertices_to_plot[:, 0].min()) * 0.5
           mid_y = (vertices_to_plot[:, 1].max() + vertices_to_plot[:, 1].min()) * 0.5
           mid_z = (vertices_to_plot[:, 2].max() + vertices_to_plot[:, 2].min()) * 0.5

           ax.quiver(mid_x, mid_y, mid_z, max_range*1.05 , 0, 0, color='red', arrow_length_ratio=0.1, label = "X")  
           ax.quiver(mid_x, mid_y, mid_z, 0, max_range*1.05, 0, color='green', arrow_length_ratio=0.1, label = "Y")  
           ax.quiver(mid_x, mid_y, mid_z, 0, 0, max_range*1.05, color='blue', arrow_length_ratio=0.1, label = "Z")  
           
           # Adjust the view for each subplot
           ax.view_init(elev=view[0], azim=view[1])
           ax.axis('off')
           plt.legend()
         
       plt.show(block=block)
                    
        
    def plot_vertices(self, transformed = False, block = False):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')
        
        vertices_to_plot = self.vertices_transformed if transformed else self.vertices
        tris = self.triangles
        
        # Create the plot with a grey surface and no lines
        #ax.plot_trisurf(vertices_to_plot[:, 0], vertices_to_plot[:, 1], vertices_to_plot[:, 2], triangles=tris, 
         #               color='grey', edgecolor='none', alpha=0.8)
        
        num_triangles_to_exclude = min(10, len(tris))
        tris_grey = tris[num_triangles_to_exclude:]  # Triangles à dessiner en gris (à partir du 11e)
        tris_blue = tris[:num_triangles_to_exclude]  # Les 10 premiers triangles à dessiner en bleu
     
        ax.plot_trisurf(vertices_to_plot[:, 0], vertices_to_plot[:, 1], vertices_to_plot[:, 2],
                        triangles=tris_grey, color='grey', edgecolor='none', alpha=0.7)
     
        ax.plot_trisurf(vertices_to_plot[:, 0], vertices_to_plot[:, 1], vertices_to_plot[:, 2],
                        triangles=tris_blue, color='blue', edgecolor='none', alpha=0.9)
        
        # Add XYZ axes for better orientation
        max_range = np.array([vertices_to_plot[:, 0].max() - vertices_to_plot[:, 0].min(), 
                              vertices_to_plot[:, 1].max() - vertices_to_plot[:, 1].min(), 
                              vertices_to_plot[:, 2].max() - vertices_to_plot[:, 2].min()]).max() / 2.0
        mid_x = (vertices_to_plot[:, 0].max() + vertices_to_plot[:, 0].min()) * 0.5
        mid_y = (vertices_to_plot[:, 1].max() + vertices_to_plot[:, 1].min()) * 0.5
        mid_z = (vertices_to_plot[:, 2].max() + vertices_to_plot[:, 2].min()) * 0.5

        ax.quiver(mid_x, mid_y, mid_z, max_range, 0, 0, color='blue', arrow_length_ratio=0.1, label = "X")  
        ax.quiver(mid_x, mid_y, mid_z, 0, max_range, 0, color='green', arrow_length_ratio=0.1, label = "Y")  
        ax.quiver(mid_x, mid_y, mid_z, 0, 0, max_range, color='red', arrow_length_ratio=0.1, label = "Z")  
    
        ax.axis('off')
        plt.legend()
        plt.show(block=block)
        
    def equivalentVolume(self, coef):
        
       scaled_vertices = self.vertices * coef

       volume = 0.0

    
       for triangle in self.triangles:
           # Récupérer les sommets du triangle
           p1 = scaled_vertices[triangle[0]]
           p2 = scaled_vertices[triangle[1]]
           p3 = scaled_vertices[triangle[2]]
           
           # Les vecteurs des sommets par rapport à l'origine sont les sommets eux-mêmes
           v0 = p1
           v1 = p2
           v2 = p3
           
           # Calcul du volume du tétraèdre (1/6 du produit mixte)
           dV = np.dot(v0, np.cross(v1, v2)) / 6.0
           volume += dV
    
       # Prendre la valeur absolue du volume total
       volume = abs(volume)
       print(f"Volume équivalent de l'astéroïde : {volume} unités cubiques")
       return volume
        
    def equivalentVolRadius(self, coef):
       eqV = self.equivalentVolume(coef)
       return (3/ (4 * np.pi) * eqV)**(1/3)
        
    def equivalentSurfaceDiameter(self, coef=1.0):
        """
        Calcule le diamètre équivalent basé sur la surface.
        1) On scale d'abord les sommets par 'coef'
        2) On calcule la surface du maillage
        3) On en déduit le diamètre d'une sphère de même surface
        
        :param coef: Facteur d'échelle pour les sommets (similaire à la méthode equivalentVolume).
        :return: Diamètre équivalent (float)
        """
        # 1) Mise à l'échelle des sommets
        scaled_vertices = self.vertices * coef
        
        # 2) Calcul de la surface totale
        total_area = 0.0
        for triangle in self.triangles:
            p1 = scaled_vertices[triangle[0]]
            p2 = scaled_vertices[triangle[1]]
            p3 = scaled_vertices[triangle[2]]
            
            # Vecteurs du triangle
            v1 = p2 - p1
            v2 = p3 - p1
            
            # Aire du triangle = 1/2 * norme( v1 x v2 )
            area_triangle = 0.5 * np.linalg.norm(np.cross(v1, v2))
            total_area += area_triangle
        
        # 3) Calcul du diamètre équivalent
        # Rayon d'une sphère de même surface : R = sqrt( total_area / (4*pi) )
        # Diamètre = 2 * R
        eq_diameter = 2.0 * np.sqrt(total_area / (4.0 * np.pi))
        
        print(f"Surface équivalente de l'astéroïde : {total_area:.4f} unités carrées")
        print(f"Diamètre équivalent (basé sur la surface) : {eq_diameter:.4f} unités")
        
        return eq_diameter
    
class Occultation:
    def __init__(self, event):
        self.event = event
        self.origineEvent = event.findFirstPositif()
        
        self.cords = {}
        self.cordsErr = {}
        self.cordsStatue = {} # 0 - positif 1 - negatif 2 - only one point
        self.cordWeight = {}
        
        dec = self.deg2rad(event.star["dec"]) #rad 
        ra = self.deg2rad(event.star["ra"]) # rad
        
        self.mean_direction = None
        self.a = None
        
        self.t_reference = None
        
        self.sEps = np.asarray([ -np.sin(dec) * np.cos(ra),
                                 -np.sin(dec) * np.sin(ra),
                                 np.cos(dec)])
        
        self.sNu = np.asarray([ np.sin(ra),
                                -np.cos(ra),
                                0])
        self.speedAvg = None
        
        self.__buildCords()
        
    def getAst(self):
        return self.event.getAst()
    
    def __len__(self):
        return len(self.event)        
        
    def __str__(self):
        
        strr = "Date "  + str(self.event.date.value) + " : \n"
        
        for c in self.cords:
            
            strr += " ID: " + str(c) + " : \n"
            strr += "  D: " + str(self.cords[c][0, 0])  + " sEps " + str(self.cords[c][0, 1]) + " sNu \n"
            strr += "  R: " + str(self.cords[c][1, 0])  + " sEps " + str(self.cords[c][1, 1]) + " sNu \n"
        
        return strr
    
    def date(self):
        return self.event.date.to_value("iso")
        
    
    def t_ref(self):
        
        if self.t_reference == None:
     
            self.t_reference = 0
            cpt = 0
            for ev in self.event:
                if ev["status"] == 0:
                    cpt += 1
                    self.t_reference += ev["d"].jd
           
            self.t_reference /= cpt 
           
            #self.t_reference = np.mean([ev["d"].jd for ev in self.event])
            
            
        return self.t_reference
    
    def avgSpeed(self):
        
        if not isinstance(self.speedAvg, type(None)):
            return self.speedAvg
        
        dv = np.zeros((3))
        
        index = self.event.index
        self.event.index = 0
        
        for ev in self.event:
            
            if ev["status"] != 2:
                x_d, x_r, x_d_err, x_r_err = self.obsPosEv(ev)
                v_d = self.speed(self.event.ast["id"], ev["d"].to_value("jd"))
                v_r = self.speed(self.event.ast["id"], ev["r"].to_value("jd")) 
                dv = dv + ( v_r + v_d ) / 2
        
        self.event.index = index
        
        self.speedAvg = dv
        
        return self.speedAvg
    
    def __buildCords(self):
        
 
        t_refer = self.t_ref() * 86400  # en secondes
  
        for ev in self.event:

            x_d, x_r, x_d_err, x_r_err = self.obsPosEv(ev)

            if ev["status"] != 2:
                v_d = self.speed(self.event.ast["id"], ev["d"].to_value("jd"))
                v_r = self.speed(self.event.ast["id"], ev["r"].to_value("jd"))
                delta_v = ( v_r + v_d ) / 2
            
                t_d = ev["d"].jd * 86400  # en secondes
                delta_t_d = t_d - t_refer  # en secondes
               
                t_r = ev["r"].jd * 86400  # en secondes
                delta_t_r = t_r - t_refer  # en secondes
                
                da = (v_d - v_r) / (delta_t_d - delta_t_r)
                
                dad = 0.5 * da * delta_t_d**2
                dar = 0.5 * da * delta_t_r**2
                
                
                x_plus_d = x_d + delta_v * delta_t_d 
                xi_d = np.dot(self.sEps, x_plus_d) / 1000
                eta_d = np.dot(self.sNu, x_plus_d) / 1000
                    
                x_plus_r = x_r + delta_v * delta_t_r
                xi_r = np.dot(self.sEps, x_plus_r) / 1000
                eta_r = np.dot(self.sNu, x_plus_r) / 1000
    
            
                if not np.isnan(da).any():
                    x_plus_d += dad
                    x_plus_r += dar
                    
            elif Event.isDmissing(ev):
                eta_d = None
                xi_d = None
                delta_v = self.speed(self.event.ast["id"], ev["r"].to_value("jd"))
                t_r = ev["r"].jd * 86400  # en secondes
                delta_t_r = t_r - t_refer  # en secondes
                x_plus_r = x_r + delta_v * delta_t_r
                xi_r = np.dot(self.sEps, x_plus_r) / 1000
                eta_r = np.dot(self.sNu, x_plus_r) / 1000
            else:
                
                eta_r = None
                xi_r = None
                delta_v = self.speed(self.event.ast["id"], ev["d"].to_value("jd"))
                t_d = ev["d"].jd * 86400  # en secondes
                delta_t_d = t_d - t_refer  # en secondes
                x_plus_d = x_d + delta_v * delta_t_d
                xi_d = np.dot(self.sEps, x_plus_d) / 1000
                eta_d = np.dot(self.sNu, x_plus_d) / 1000
                
            
            # ========== Calcul des erreurs ===============
            # Temps avec erreurs
            if not Event.isDmissing(ev):
                t_d_err_p = (ev["d"] + ev["d_err"]).jd * 86400
                t_d_err_n = (ev["d"] - ev["d_err"]).jd * 86400
                
                delta_t_d_err_p = t_d_err_p - t_refer
                delta_t_d_err_n = t_d_err_n - t_refer
                
                x_plus_d_err_p = x_d + delta_v * delta_t_d_err_p
                x_plus_d_err_n = x_d + delta_v * delta_t_d_err_n
                
                if not np.isnan(da).any():
                    dad_err_p = 0.5 * da * delta_t_d_err_p**2
                    dad_err_n = 0.5 * da * delta_t_d_err_n**2
                    
                    x_plus_d_err_p += dad_err_p
                    x_plus_d_err_n += dad_err_n
                
                xi_d_err_p = np.dot(self.sEps, x_plus_d_err_p) / 1000
                eta_d_err_p = np.dot(self.sNu, x_plus_d_err_p) / 1000
                xi_d_err_n = np.dot(self.sEps, x_plus_d_err_n) / 1000
                eta_d_err_n = np.dot(self.sNu, x_plus_d_err_n) / 1000
           
            else:
                xi_d_err_p = None
                eta_d_err_p = None
                xi_d_err_n = None
                eta_d_err_n = None
        
            if not Event.isRmissing(ev):
                t_r_err_p = (ev["r"] + ev["r_err"]).jd * 86400
                t_r_err_n = (ev["r"] - ev["r_err"]).jd * 86400
            
                
                delta_t_r_err_p = t_r_err_p - t_refer
                delta_t_r_err_n = t_r_err_n - t_refer
                
               
                x_plus_r_err_p = x_r + delta_v * delta_t_r_err_p
                x_plus_r_err_n = x_r + delta_v * delta_t_r_err_n
                
                if not np.isnan(da).any():
                    dar_err_p = 0.5 * da * delta_t_r_err_p**2
                    dar_err_n = 0.5 * da * delta_t_r_err_n**2
                
                    x_plus_r_err_p += dar_err_p
                    x_plus_r_err_n += dar_err_n
            
        
                xi_r_err_p = np.dot(self.sEps, x_plus_r_err_p) / 1000
                eta_r_err_p = np.dot(self.sNu, x_plus_r_err_p) / 1000
                xi_r_err_n = np.dot(self.sEps, x_plus_r_err_n) / 1000
                eta_r_err_n = np.dot(self.sNu, x_plus_r_err_n) / 1000
           
            else:
                xi_r_err_p = None
                eta_r_err_p = None
                xi_r_err_n = None
                eta_r_err_n = None
            
            # Stockage des erreurs projetées (en kilomètres)

            # Stockage des cordes (en kilomètres)
            self.cords[self.event.getId(ev)] = np.array([[xi_d, eta_d],
                                                         [xi_r, eta_r]])

           
            self.cordsErr[self.event.getId(ev)] = {
                 'd_err_p': np.array([xi_d_err_p, eta_d_err_p]),
                 'd_err_n': np.array([xi_d_err_n, eta_d_err_n]),
                 'r_err_p': np.array([xi_r_err_p, eta_r_err_p]),
                 'r_err_n': np.array([xi_r_err_n, eta_r_err_n])
             }
            
            self.cordsStatue[self.event.getId(ev)] = ev["status"]
            
            # if (self.cords[self.event.getId(ev)][0, : ] == self.cords[self.event.getId(ev)][1, : ]).all():
            #     self.cordsStatue[self.event.getId(ev)] = 1
            # else:
            #     self.cordsStatue[self.event.getId(ev)] = 0
            
            self.cordWeight[self.event.getId(ev)] = 1.0
                
        self.computeTendances()
                
            
    def computeTendances(self):
        
        positive_vectors = []
        self.a = 0
        for cordId in self.cords:
            status = self.cordsStatue[cordId]
            if status == 0:
                cord = self.cords[cordId]
                vector = cord[1] - cord[0]
                
                self.a += vector[0] / vector[1]
                
                norm = np.linalg.norm(vector)
                if norm != 0:
                    unit_vector = vector / norm
                    positive_vectors.append(unit_vector)
    
        if len(positive_vectors) == 0:
            print("no positive cords to compute the tendance")
            return
        self.a /= len(positive_vectors)
        # Direction moyenne des cordes positives
        self.mean_direction = np.mean(positive_vectors, axis=0)
        self.mean_direction /= np.linalg.norm(self.mean_direction)  # Normalisation
                
    def getAsteroidPos(self, astId, time):
        jdTime = time.to_value('jd')
        obj = Horizons(id=astId, id_type='smallbody', location='500@399', epochs=jdTime)
        vec = obj.vectors(refplane='earth')  # Coordonnées géocentriques équatoriales
    
        x = vec['x'][0] * u.au.to(u.m)  # Conversion en mètres
        y = vec['y'][0] * u.au.to(u.m)
        z = vec['z'][0] * u.au.to(u.m)
        return np.array([x, y, z])
    
    def animateFondPlan(self, invertXaxis = True):
        
        fig, ax = plt.subplots()
        
        
            
        v = self.speed(self.event.ast["id"], self.origineEvent["d"].to_value("jd"))
        v_projected = self.posProjected(v)
        
        
        def animation(i):
            
            time = i / 100 #s 
        
            ax.clear()
        
       
            lat = self.hAngle2Degree(self.origineEvent["lat"]) # degree
            long = self.hAngle2Degree(self.origineEvent["long"]) # degree
            
            
            origine = self.obsPos(lat, long, self.origineEvent["alt"], self.origineEvent["d"])
            origineProj = self.posProjected(origine)
            
            
            for ev in self.event:
                
                lat = self.hAngle2Degree(ev["lat"]) # degree
                long = self.hAngle2Degree(ev["long"]) # degree
                
                pos = self.obsPos(lat, long, ev["alt"], self.origineEvent["d"])
                pos = self.posProjected(pos)
        
                ax.scatter(pos[1], pos[0])
                ax.text(pos[1], pos[0], str(self.event.getId(ev)))
    
            ast_pos = origineProj + v_projected * time
            
            ax.scatter(ast_pos[1], ast_pos[0])
            ax.xaxis.set_inverted(invertXaxis)
            fig.canvas.draw()
            
            plt.pause(0.001)
     
        animation(200)
        
        
        #for i in range(8000):
         #   animation(i)
        
       # anim = FuncAnimation(fig, animation, frames = None)

       
    def correctPosShift(self, x_projected, speed_projected):
        
        
        origine = self.getOriginePos()
        origine = self.posProjected(origine)
        
        vec = x_projected - origine
        norm_vec = np.linalg.norm(vec)
        
        eps = speed_projected[1] * norm_vec / (speed_projected[1] - speed_projected[0])
        nu = np.sqrt(norm_vec**2 - eps**2)
        
        return np.asarray([eps, nu])
        
    
    

    def correctEvFromTimeShift(self, ev):
        
        if ev == self.origineEvent:
            return
        
        lat = self.hAngle2Degree(ev["lat"]) # degree
        long = self.hAngle2Degree(ev["long"]) # degree
        
        
        
        x_obs = self.obsPos(lat, long, ev["alt"], self.origineEvent["d"])
        x_origine = self.getOriginePos()
        
        x_obs_projected = self.posProjected(x_obs)
        x_origine_projected = self.posProjected(x_origine)
        
        v = self.speed(self.event.ast["id"], self.origineEvent["d"].to_value("jd"))
        v_projected = self.posProjected(v)
            
        dist = (x_obs_projected - x_origine_projected)
        v_projected_along_dist = np.dot(dist, v_projected) / np.linalg.norm(dist)
        
        dt = np.linalg.norm(dist) / v_projected_along_dist
        
    
        newDtime = self.time2decimal(ev["d"].value) - dt / 3600
        newRtime = self.time2decimal(ev["r"].value) - dt / 3600
       
        #-----------------
        
        if newDtime < 0:

           date = ev["d"].value.split(" ")[0].split("-")
           date[-1]=str(int( ev["d"].value.split(" ")[0].split("-")[-1]) - 1)
           
           newDtime = '-'.join(date) + " " + self.decimal2time(24 + newDtime)
           
        else:   
            newDtime = ev["d"].value.split(" ")[0] + " " + self.decimal2time(newDtime)
            
        if newRtime < 0:
            date = ev["r"].value.split(" ")[0].split("-")
            date[-1]=str(int( ev["r"].value.split(" ")[0].split("-")[-1]) - 1)
               
            newRtime = '-'.join(date) + " " + self.decimal2time(24 + newRtime)
        else:   
            newRtime = ev["r"].value.split(" ")[0] + " " + self.decimal2time(newRtime)
       
        #-----------------

        ev["d"] = Time(newDtime, format = "iso")
        ev["r"] = Time(newRtime, format = "iso")
        
        
    
    def getOriginePos(self):
        olat = self.hAngle2Degree(self.origineEvent["lat"]) # degree
        olong = self.hAngle2Degree(self.origineEvent["long"]) # degree
        
        return self.obsPos(olat, olong, self.origineEvent['alt'], self.origineEvent["d"])
    
    def posProjectedEv(self, ev):
        x_d, x_r, _, _ = self.obsPosEv(ev)
        return self.posProjected(x_d),  self.posProjected(x_r)
    
    def posProjected(self, x):
        return np.asarray([np.dot(self.sEps, x), np.dot(self.sNu, x)])
        
    def avgTime(self):
        
        
        delta = 0
        cpt = 0
        
        index = self.event.index
        self.event.index = 0                    
                    
        for ev in self.event:
            if ev["isPositif"]:
                delta += (self.time2decimal(ev["r"].value) - self.time2decimal(ev["d"].value)) * 3600
                cpt += 1
                
        self.event.index = index
        
        return delta / cpt
    
    def obsPosEv(self, ev):
        lat = self.hAngle2Degree(ev["lat"]) # degree
        long = self.hAngle2Degree(ev["long"]) # degree
        
        if not Event.isDmissing(ev):
            x_d = self.obsPos(lat, long, ev["alt"], ev["d"])
            x_d_err_m = self.obsPos(lat, long, ev["alt"], ev["d"] - ev['d_err'])
            x_d_err_p = self.obsPos(lat, long, ev["alt"], ev["d"] + ev['d_err']) 
        else:
            x_d = None
            x_d_err_m = None
            x_d_err_p = None
        
        if not Event.isRmissing(ev):
            x_r = self.obsPos(lat, long, ev["alt"], ev["r"])
            x_r_err_m = self.obsPos(lat, long, ev["alt"], ev["r"] - ev['r_err'])
            x_r_err_p = self.obsPos(lat, long, ev["alt"], ev["r"] + ev['r_err']) 
        else:
            x_r = None
            x_r_err_m = None
            x_r_err_p = None
        
        return x_d, x_r, [x_d_err_m, x_d_err_p], [x_r_err_m, x_r_err_p]
    

        
    def obsPos(self, lat, long, height, time):
        
        time = time.value    
    
        R = EarthLocation(lat=lat*u.deg, lon=long * u.deg, height=height*u.m).geocentric  # m
        R = np.asarray([R[0].value, R[1].value, R[2].value]) # m
        R = np.sqrt(np.sum(R*R)) # m distance to the earth center
        
        theta = self.LST(time, lat, long )

        
        theta = theta.rad
        
        phi = self.deg2rad(lat)
        
        x = np.asarray([ R * np.cos(phi) * np.cos(theta),
                         R * np.cos(phi) * np.sin(theta),
                         R * np.sin(phi)])
        
        return x
    
    
    def obsPos2(self, lat, long, height, time):
        long = long
        location = EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=height*u.m)
        gcrs = location.get_gcrs(obstime=time)
        x = gcrs.cartesian.x.value  # en mètres
        y = gcrs.cartesian.y.value
        z = gcrs.cartesian.z.value
        return np.array([x, y, z])
    
    
    def speed(self, astId, jdTime):
        
        
        h = Horizons(id=astId, id_type = "smallbody", location= "geocentric", epochs = jdTime)
       
        v = np.asarray([  self.aud2ms(h.vectors(refplane='earth')["vx"].data[0]),
                            self.aud2ms(h.vectors(refplane='earth')["vy"].data[0]),
                            self.aud2ms(h.vectors(refplane='earth')["vz"].data[0] )]) # m / s
        
        return -v
    
    
    def aud2ms(self, aud):
        return (aud * u.AU / u.day).to(u.m / u.s).value
    
    def deg2rad(self, deg):
        return deg * np.pi / 180
    
    def rad2deg(self, rad):
        return rad * 180 / np.pi
    
    def hAngle2Degree(self, hAngle):
        
        
        
        s = [x for x in hAngle.split(" ") if x != ""]
        
        if len(s) == 4 and (s[0] == "+" or s[0] == "-"):
            s[1] = s[0] + s[1]
            s = s[1:]
        
        sign = s[0][0]
        
        if sign == "-":
            return float(s[0]) - float(s[1])/60 - float(s[2]) / 3600 
        else:
            return float(s[0]) + float(s[1])/60 + float(s[2]) / 3600 
        
   
        
    def hAngle2Rad(self, hAngle):
        return self.deg2rad( self.hAngle2Degree(hAngle) )
    
    def time2decimal(self, date):
        hms = date.split(" ")[-1].split(':')
        
   
        h = float(hms[0])
        m = float(hms[1]) / 60
        s = float(hms[2]) / 3600
        
        return h + m + s

    def decimal2time(self, time):
        
     
        h = int(time)
        minute = int((time - h) * 60)
        second = ((time - h) * 60 - minute) * 60
        
        return str(h) + ":" + str(minute) + ":" + str(second)
   
 

    def LST(self, datetime, latDecimal, longDecimal):
       
      
        
        observing_location = EarthLocation(lat=latDecimal*u.deg, lon=longDecimal*u.deg)
        observing_time = Time(datetime, format = "iso", scale='ut1')
        
        temps_sideral_local = observing_time.sidereal_time('apparent', longitude=Longitude(longDecimal, unit='deg'))
        
        return temps_sideral_local
    
    def plot(self, invertXaxis=False):
        fig, ax = plt.subplots()
        
        for c in self.cords:
            ax.plot(self.cords[c][:, 1], self.cords[c][:, 0], marker='o', linestyle='-', label=f'Observer {c}')
            ax.text(self.cords[c][0, 1], self.cords[c][0, 0], str(c))
    
        ax.set_xlabel('Coordonnée η (km)')
        ax.set_ylabel('Coordonnée ξ (km)')
        ax.set_title('Projection des cordes dans le plan fondamental (correction du décalage)')
        ax.grid(True)
        ax.legend()
        if invertXaxis:
            ax.invert_xaxis()
        plt.show(block=False)
    
    def plot3dChart(self, shouldDrawSphere = True, shouldDrawPlan = True):
        
        
        def drawSphere(ax):
            
            t_avg = np.mean([self.time2decimal(ev["d"].value) for ev in self.event if ev["isPositif"]])
     
            zro = self.obsPos(0, 0, 0, Time(self.event.date.value.split(' ')[0] + " " + self.decimal2time( t_avg ), format = "iso"))
            maxx = np.linalg.norm(zro)
            zro = zro/maxx
            
            bm =  Image.open(file_path + 'Equirectangular_projection_SW.jpg')

            bm = np.array(bm.resize([int(d/5) for d in bm.size]))/256.
            
            #bm = np.fliplr(bm)

            # coordinates of the image - don't know if this is entirely accurate, but probably close
            lons = np.linspace(0, 360, bm.shape[1]) * np.pi/180 
            lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 
            
            # repeat code from one of the examples linked to in the question, except for specifying facecolors:
            x = np.outer(np.cos(lons), np.cos(lats)).T
            y = np.outer(np.sin(lons), np.cos(lats)).T
            z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
            
   
            theta_rotation = np.pi + np.arccos(zro[0]) * np.sign(np.arcsin(zro[1]))

            x_rot = x * np.cos(theta_rotation) - y * np.sin(theta_rotation)
            y_rot = x * np.sin(theta_rotation) + y * np.cos(theta_rotation)

            
            ax.plot_surface(x_rot, y_rot, z, facecolors = bm, edgecolor='none', alpha= 0.5)
            ax.set_aspect('auto') 
            

        def drawPlan(ax, point1, point2, point3):
            
            v1 = point2 - point1
            v2 = point3 - point1
            
            normal = np.cross(v1, v2)
            
            if np.linalg.norm(normal) == 0:
                print("Les trois points sont colinéaires et ne définissent pas un plan unique.")
            else:
                # Équation du plan : ax + by + cz + d = 0
                a, b, c = normal
                d = -np.dot(normal, point1)
            
                xx, yy = np.meshgrid(range(-1, 3), range(-1, 3))
            
                if c != 0:
                    zz = (-a * xx - b * yy - d) / c
                else:
                    # Si c est nul, le plan est vertical, on calcule x en fonction de y et z
                    zz = np.linspace(-1, 3, 4)
                    yy, zz = np.meshgrid(range(-1, 3), zz)
                    xx = (-b * yy - c * zz - d) / a
                
            ax.plot_surface(xx, yy, zz, alpha=0.5)
        
        t_avg = np.mean([self.time2decimal(ev["d"].value) for ev in self.event if ev["isPositif"]])
        zro = self.obsPos(0, 0, 0, Time(self.event.date.value.split(' ')[0] + " " + self.decimal2time( t_avg ), format = "iso"))
        maxx = np.linalg.norm(zro)
        
        
        
        ra = self.deg2rad(102.95)
        dec = self.deg2rad(-34.39)
        x_star_az = np.cos(dec) * np.cos(ra)
        y_star_az = np.cos(dec) * np.sin(ra)
        z_star_az = np.sin(dec)
        
        ra = np.radians(self.event.star["ra"])
        dec = np.radians(self.event.star['dec'])
        
        
       
        v = self.speed(self.event.ast["id"], self.origineEvent['d'].to_value("jd"))
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if shouldDrawSphere:
            drawSphere(ax)
        ax.scatter(zro[0] / maxx, zro[1] / maxx, zro[2] / maxx, color = "orange")
        
        ax.quiver(0, 0, 0, 1, 0, 0, color='b', label='X (Point Vernal)')  # X axis (vers point vernal)
        ax.quiver(0, 0, 0, 0, 1, 0, color='g', label='Y')  # Y axis
        ax.quiver(0, 0, 0, 0, 0, 1, color='k', label='Z (Pôle Nord)')  # Z axis (projection du pôle nord)
        
        x_star = np.cos(dec) * np.cos(ra)
        y_star = np.cos(dec) * np.sin(ra)
        z_star = np.sin(dec)
        ax.scatter(x_star, y_star, z_star, color='y', s=100, label='Étoile')
        ax.plot([0, x_star], [0, y_star], [0, z_star], color='m', linestyle='--')

        
        # Projection de l'étoile sur le plan XY
        x_proj = np.cos(dec) * np.cos(ra)
        y_proj = np.cos(dec) * np.sin(ra)
        z_proj = 0
        
        # Tracer l'ascension droite (RA) comme un arc de cercle entre X et la projection sur le plan XY
        ra_angle_deg = np.degrees(ra)
        arc_ra = Arc((0, 0), 0.5, 0.5, theta1=0, theta2=ra_angle_deg, color='c', linestyle='--')
        ax.add_patch(arc_ra)
        art3d.pathpatch_2d_to_3d(arc_ra, z=0, zdir='z')
        
        # Tracer la déclinaison (Dec) comme un arc de cercle entre l'axe allant du centre du repère à l'étoile et le plan XY
        dec_radius = np.sqrt(x_proj**2 + y_proj**2)  # Rayon de la projection sur le plan XY
        dec_angle_deg = np.degrees(dec)
        
        if dec_angle_deg > 0:
            arc_dec = Arc((0, 0, 0), 1, 1, theta1=0, theta2=dec_angle_deg, color='m', linestyle='--')
        else:
            arc_dec = Arc((0, 0, 0), 1, 1, theta1=dec_angle_deg, theta2=0, color='m', linestyle='--')
        
        ax.add_patch(arc_dec)
        art3d.pathpatch_2d_to_3d(arc_dec, z=0, zdir='y')
        
        
        arc_equateur = Arc((0, 0, 0), 2, 2, theta1=2*np.pi, theta2=0, color='c', linestyle='--')
        ax.add_patch(arc_equateur)
        art3d.pathpatch_2d_to_3d(arc_equateur, z=0, zdir='z')

        eps_x = self.sEps[0]
        eps_y = self.sEps[1]
        eps_z = self.sEps[2]
        ax.quiver(0, 0, 0, eps_x, eps_y, eps_z, color='r', label='Vecteur eps')
        
        nu_x = self.sNu[0]
        nu_y = self.sNu[1]
        nu_z = self.sNu[2]
        ax.quiver(0, 0, 0, nu_x, nu_y, nu_z, color='grey', linestyle='-', label='Vecteur nu')
        
        if shouldDrawPlan:
            drawPlan(ax, np.array([eps_x, eps_y, eps_z]), np.array([nu_x, nu_y, nu_z]), np.array([0,0,0]))
        
        # Représenter le vecteur v dans le plan eps-nu
        if v is not None:
            # Calculer la composante de v dans le plan défini par eps et nu
            v_proj = self.posProjected(v)
            
            maxx = np.max(np.abs(v))
            # Tracer le vecteur v_total
          
            ax.quiver(0, 0, 0, v[0] / maxx, v[1] / maxx, v[2] / maxx, color='orange', label='Vecteur v')
          
            v_eps = v_proj[0] * np.array([eps_x, eps_y, eps_z])
            v_nu = v_proj[1] * np.array([nu_x, nu_y, nu_z])
            v_total = v_eps + v_nu
            
            maxx = np.max(np.abs(v_total))
            ax.quiver(0, 0, 0, v_total[0] / maxx, v_total[1]/ maxx, v_total[2]/ maxx, color='purple', label='Vecteur v (dans le plan eps-nu)')
        
        for ev in self.event:
            
            x_d,_ , _, _= self.obsPosEv(ev)
            x_d_proj = self.posProjected(x_d)
            
            x_d_eps = x_d_proj[0] * np.array([eps_x, eps_y, eps_z])
            x_d_nu = x_d_proj[1] * np.array([nu_x, nu_y, nu_z])
            x_total = x_d_eps + x_d_nu
            
            maxx = np.linalg.norm(x_d)
            ax.scatter(x_d[0] / maxx, x_d[1] / maxx, x_d[2] / maxx, color = "black")
            maxx = np.linalg.norm(x_total)
            ax.scatter(x_total[0] / maxx, x_total[1] / maxx, x_total[2] / maxx, color = "blue")

            
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        

        ax.legend()
        
     
        plt.show(block=False)
        
        return fig, ax
    
class Polygon:
    
    def __init__(self, vertices, edges):
       
        self.segments = []
         
        for edge in edges:
             p1, p2 = edge
             self.segments.append( np.asarray([ vertices[p1], vertices[p2] ] ) )

    def __point2line(self, point1, point2):
        y1, x1 = point1
        y2, x2 = point2
        
        if x1 == x2:  # Cas de droite verticale
            return None, x1  # Utilise 'None' pour indiquer une pente infinie
        
        a = (y2 - y1) / (x2 - x1)
        b = y1 - a * x1
        return a, b

    def __lineIntersectLine(self, a1, b1, a2, b2):
        if a1 == a2:  # Droites parallèles
            return None
        if a1 is None:  # La première droite est verticale
            x = b1
            y = a2 * x + b2
        elif a2 is None:  # La deuxième droite est verticale
            x = b2
            y = a1 * x + b1
        else:
            x = (b2 - b1) / (a1 - a2)
            y = a1 * x + b1
        return np.asarray([x, y])
    
    def __isInterInSeg(self, segment, point):
        y1, x1 = segment[0]
        y2, x2 = segment[1]
        x, y = point
        
        in_bounds_x = min(x1, x2) <= x <= max(x1, x2)
        in_bounds_y = min(y1, y2) <= y <= max(y1, y2)
        
        #print(y1, x1, y2, x2, y, x, in_bounds_x and in_bounds_y)
        
        return in_bounds_x and in_bounds_y 

    def intersect(self, point1, point2):
        
        intersection = []
        
        a, b = self.__point2line(point1, point2)
        #print(point1, point2, self.segments)

        for seg in self.segments:
            
            aseg, bseg = self.__point2line(seg[0, :], seg[1, :])
            inter = self.__lineIntersectLine(a, b, aseg, bseg)
            
                        
            #print(seg[0,:], seg[1,:], point1, point2)
            
            if not isinstance(type(inter), type(None)) and self.__isInterInSeg(seg, inter):
                intersection.append(inter)
        
        return intersection

class OcFitInputParser:
    def __init__(self, path):
        
        self.nSolution = 0
        self.lamb = []
        self.beta = []
        self.p = []
        self.zeroT = []
        self.initPhase = []
        self.p_err_p = []
        self.p_err_m = []
        self.l_err_p = []
        self.l_err_m = []
        self.b_err_p = []
        self.b_err_m = []
        
        self.l_opt = []
        self.b_opt = []
        self.p_opt = []
        
        self.cordToDelete = {}
        self.cordToShift = {}
        self.weight = {}
        self.optimizer = {}
        self.multi = {}
        self.centroid = {}
        self.size = {}
        
        self.__parse(path)
        
    def __removeEmptyFromList(self, li):
        return [x for x in li if len(x) > 0]
    
    def __parseSolution(self, solution):
        
        for i in range(self.nSolution):
            
            s = solution[1 + i].split('\n')
            
            for i, el in enumerate(s): #remove comments
                s[i] = el.split("#")[0]
            
            self.lamb.append(float(s[0].split("+")[0]))
            self.l_err_p.append( float( s[0].split("+")[1].split("-")[0]) )
            self.l_err_m.append( float( self.__removeEmptyFromList( s[0].split("-")[1].replace(" ", "\t").split("\t") )[0] ) )
            self.l_opt.append( float(self.__removeEmptyFromList( s[0].split("-")[-1].replace(" ", "\t").split("\t") )[1] ) )
            
            self.beta.append(float(s[1].split("+")[0]))
            self.b_err_p.append( float( s[1].split("+")[1].split("-")[0]) )
            self.b_err_m.append( float( self.__removeEmptyFromList( s[1].split("-")[1].replace(" ", "\t").split("\t") )[0] ) )
            self.b_opt.append( float( self.__removeEmptyFromList( s[1].split("-")[-1].replace(" ", "\t").split("\t") )[1]  ) )
            
            self.p.append(float(s[2].split("+")[0]))
            self.p_err_p.append( float( s[2].split("+")[1].split("-")[0]) )
            self.p_err_m.append( float( self.__removeEmptyFromList( s[2].split("-")[1].replace(" ", "\t").split("\t") )[0] ) )
            self.p_opt.append( float(self.__removeEmptyFromList( s[2].split("-")[-1].replace(" ", "\t").split("\t") )[1] ) )
            
            self.zeroT.append( float(self.__removeEmptyFromList(s[3].replace(" ", "\t").split("\t"))[0]))
            self.initPhase.append( float(self.__removeEmptyFromList(s[3].replace(" ", "\t").split("\t"))[1]) )
            
            
    def __parseOneColumn(self, inp):
        lines = inp.split("\n")[1:]
        
        elements = []
        for line in lines:
            
            if len(line) == 0:
                continue
            el = self.__removeEmptyFromList(line.replace('\t', " ").split(" "))
            elements += el
        
        return elements
            
    def __parseTwoColumn(self, inp):
        lines = inp.split("\n")[1:]
        
        
        keys = []
        els = []
    
        for line in lines:
            
            if len(line) == 0:
                continue
            
            print(line)
            
            elements = self.__removeEmptyFromList(line.replace('\t', " ").split(" "))
            keys.append(elements[0]) 
            els.append(elements[1])
            
        return keys, els
            
    def __parseThreeColumns(self, inp):
        lines = [f for f in inp.split("\n") if f != ""]

        keys = []
        els1 = []
        els2 = []


        for line in lines:
  
            if len(line.replace('\t','').replace(" ", "")) == 0:
                continue
            
            elements = self.__removeEmptyFromList(line.replace('\t', " ").split(" "))

            keys.append(elements[0])
            els1.append(elements[1])
            els2.append(elements[2])
            
        return keys, els1, els2
        
    def __parseDeleteChords(self, inp):
        
        keys, els = self.__parseTwoColumn(inp)

        
        for i, key in enumerate(keys):
            if Time(key, format = "iso").to_value('iso').split(' ')[0] not in self.cordToDelete:
                self.cordToDelete[Time(key, format = "iso").to_value('iso').split(' ')[0]] = [int(els[i])]
            else:
                self.cordToDelete[Time(key, format = "iso").to_value('iso').split(' ')[0]].append(int(els[i]))
        
    def __parseShiftChords(self, inp):
        
        keys, els1, els2 = self.__parseThreeColumns(inp)
        
        for i, key in enumerate(keys):
                     
            if key not in self.cordToShift:
                self.cordToShift[key] = [[int(els1[i]), float(els2[i])]]
            else:
                self.cordToShift[key].append([int(els1[i]), float(els2[i])])
    
    def __parseWeightChords(self, inp):
    
        keys, els1, els2 = self.__parseThreeColumns(inp)
  
        for i, key in enumerate(keys):
            
            if  float(els2[i]) < 0:
                print("Error! weight for", key, "is not valide. Should be higher or equal than 0 and usually between 0 and 1 but can be higher than 1 if needed")
                sys.exit()
            
            if key not in self.weight:
                self.weight[key] = [[int(els1[i]), float(els2[i])]]
            else:
                self.weight[key].append( [int(els1[i]), float(els2[i])] )
            
    def __parseOptimizer(self, inp):
        options = inp.split("\n")
        
        options = [o.split(" ") for o in options]
        
        all_options_availbale = ["n_sub_divisions", "sub_division_max_stop"]
        
        for i in range(len(options)):
            options[i] = [o for o in options[i] if o != ""]
            
        for option in options:
            if option[0] in all_options_availbale:
                self.optimizer[option[0]] = option[1]
            
        
    def __parseMulti(self, inp):
        keys, els = self.__parseTwoColumn(inp)
        
        for i, key in enumerate(keys):
            if Time(key, format = "iso").to_value('iso').split(' ')[0] not in self.multi:
                self.multi[Time(key, format = "iso").to_value('iso').split(' ')[0]] = [int(els[i])]
            else:
                self.multi[Time(key, format = "iso").to_value('iso').split(' ')[0]].append(int(els[i]))
                
    def __parseCentroid(self, inp):
        
        keys, els1, els2 = self.__parseThreeColumns(inp) #x y

        for i, key in enumerate(keys):
            self.centroid[key] = [float(els1[i]), float(els2[i])]
        
    def __parseSize(self, inp):
        keys, els1, els2 = self.__parseThreeColumns(inp) # solution, size
        
        for i, key in enumerate(keys):
            
            self.size[key] == [None, None]
            
            solution = int(els1[i])
            val = float(els2[i])
            
            if solution == 1:
                self.size[key][0] = val
            elif solution == 2:
                self.size[key][1] = val
            elif solution == 3:
                self.size[key] = [val, val]
    
    def __parse(self, path):
        
        f = open(path, "r")
        
        inputFile = f.read()
        
        lines = inputFile.split('\n')
        lines = [l.split("#")[0] for l in lines]
        self.nSolution = int([l for l in lines[0].split(' ') if l != ""][0])
        lines = "\n".join(lines)
        
        sections = re.split(r'(?:\n\s*-{3,}\s*\n)', lines)
        self.__parseSolution(sections)
        
        
        listOfSectionPresent = []
        if 'DELETE:' in lines:
            listOfSectionPresent.append("DELETE:")
        if "SHIFT:" in lines:
            listOfSectionPresent.append("SHIFT:")
        if "WEIGHT:" in lines:
            listOfSectionPresent.append("WEIGHT:")
        if "OPTIMIZER:" in lines:
            listOfSectionPresent.append("OPTIMIZER:")
        if "MULTI:" in lines:
            listOfSectionPresent.append("MULTI:")
        if "CENTROID:" in lines:
            listOfSectionPresent.append("CENTROID:")
        if "SIZE:" in lines:
            listOfSectionPresent.append("SIZE:")
            
        
        sections = {}
        for secName in listOfSectionPresent:
            
            sections[secName] = lines.split(secName)[-1]
        
            for secName2 in listOfSectionPresent:
                
                
                if secName2 == secName:
                    continue
              
                if secName2 in sections[secName]:
                    sections[secName] = sections[secName].split(secName2)[0]
        
        if 'DELETE:' in lines:
            self.__parseDeleteChords(sections["DELETE:"])
        if "SHIFT:" in lines:
            self.__parseShiftChords(sections["SHIFT:"])
        if "WEIGHT:" in lines:
            self.__parseWeightChords(sections["WEIGHT:"])
        if "OPTIMIZER:" in lines:
            self.__parseOptimizer(sections["OPTIMIZER:"])
        if "MULTI:" in lines:
            self.__parseMulti(sections["MULTI:"])
        if "CENTROID:" in lines:
            self.__parseCentroid(sections["CENTROID:"])
        if "SIZE:" in lines:
            self.__parseSize(sections["SIZE:"])
        
        
    def getLamb(self, solution):
        
        return self.lamb[solution]
    
    def getBeta(self, solution):
        
        return self.beta[solution]
    
    def getPeriod(self, solution):
        
        return self.p[solution]
    
    def getLambErr(self, solution):
        
        return [self.l_err_m[solution], self.l_err_p[solution]]
    
    def getBetaErr(self, solution):
        
        return [self.b_err_m[solution], self.b_err_p[solution]]
    
    def getPeriodErr(self, solution):
        
        return [self.p_err_m[solution], self.p_err_p[solution]]
        
    def getInitialTime(self, solution):
        
        return self.zeroT[solution]
    
    def getInitialPhase(self, solution):
        
        return self.initPhase[solution]

class MultiParamData:
    def __init__(self, listOflistOfIdx, listOfNormaliseVector):
        
        if len(listOflistOfIdx) != len(listOfNormaliseVector):
            print("error! list of index and normalised vector should be the same lengh.", len(listOflistOfIdx), "!=" , len(listOfNormaliseVector))
            sys.exit()
        
        self.listOfListOfIdx = listOflistOfIdx
        self.listOfNormaliseVector = listOfNormaliseVector
        self.current = 0
        
    def __len__(self):
        return len(self.listOfListOfIdx)
    
    def __iter__(self):
        return self

    def __next__(self): # Python 2: def next(self)
        self.current += 1
        if self.current < len( self.listOfListOfIdx):
            return self.listOfListOfIdx[self.current]
        raise StopIteration

    def where(self, idx):
        
        for i, li in enumerate(self.listOfListOfIdx):
            if idx in li:
                return i
        
        return -1

    def keys(self):
        return self.listOfListOfIdx
    
    def vals(self):
        return self.listOfNormaliseVector
    
    def getList(self):
        
        finalLi = self.listOfListOfIdx[0]
      
        for i in range(len(self.listOfListOfIdx) -1):
            finalLi += self.listOfListOfIdx[i+1]
            
        return finalLi
    
    def toRm(self):
        
        liToRm = self.listOfListOfIdx[0][1:]
        
        for i in range(len(self.listOfListOfIdx) -1):
            liToRm += self.listOfListOfIdx[i+1][1:]
        
        return liToRm
        
        
class OccultationFit:
    def __init__(self, Occultations, model_solution1 = None, model_solution2 = None, shouldPorject = False):
        """
        Occulations is a list of object of type Occultation
        model is an object of type ModelParser
        s1 and s2 represent dictionnary of pole solutions. they should have the following keys:
            - lamb (lambda in degree)
            - err_l_m (uncertainties on lambda negatif)
            - err_l_p (uncertainties on lambda positif)
            - beta 
            - err_b_m (uncertainties on beta negatif)
            - err_b_p (uncertainties on beta positif)
        
        """
        
        self.occults = Occultations
        self.model1 = model_solution1
        self.model2 = model_solution2
        
        self.vertices1 = None
        self.vertices2 = None
        
        self.projected_vertices1 = None
        self.projected_vertices2 = None
        
        self.edge1 = None
        self.edge2 = None
        
        self.verticesMem = [None]*len(self.occults)
        self.projectedVertMem = [None]*len(self.occults)
        self.edgesMem = [None]*len(self.occults)
        
        self.param1 = None
        self.param2 = None
        
        self.bestParam1 = None
        self.bestParam2 = None
        
        self.shouldPorject = shouldPorject
    
        self.shouldMultiInit = False
        self.multi = None
        self.initialCentroids = None
        self.initialSizes = None

        self.ignorUncer = False
    

    def __initMem(self, model):
        
        for i, occ in enumerate(self.occults):
           
            t = occ.t_ref()
            lamb = model.solution["lambda"]
            beta = model.solution["beta"]
            period = model.solution["period"]
                
        
            self.__setAngle(t, model)
            self.verticesMem[i], self.projectedVertMem[i], self.edgesMem[i] = self.compute_projection2(lamb, beta, t, period, model, occ)
            
        
    
    def __spherical2cartesian(self, alpha, delta):
        
        return np.asarray([ np.cos(alpha) * np.cos(delta),
                            np.cos(delta) * np.sin(alpha),
                            np.sin(delta)])
    
    def compute_aspect_angle(self, model, occ):
        """
        Calcule l'angle d'aspect entre l'axe de rotation de l'astéroïde et la ligne de visée.
        
        Parameters:
        - model: ModelParser
            Le modèle de l'astéroïde pour lequel calculer l'angle d'aspect.
        
        Returns:
        - aspect_angle_deg: float
            L'angle d'aspect en degrés.
        """
        # 1. Calculer l'axe de rotation en coordonnées équatoriales
        lamb = model.solution['lambda']
        beta = model.solution['beta']
        alpha, delta = self.__compute_geocentric_projection(lamb, beta, occ.event.date)  # λ et β en degrés
        rotation_axis = self.__spherical2cartesian(alpha, delta)  # α et δ en radians
        rotation_axis /= np.linalg.norm(rotation_axis)
        
        # 2. Récupérer la ligne de visée en coordonnées équatoriales
        ra_star = np.deg2rad(occ.event.star["ra"])
        dec_star = np.deg2rad(occ.event.star["dec"])
        line_of_sight = self.__spherical2cartesian(ra_star, dec_star)
        line_of_sight /= -1*np.linalg.norm(line_of_sight) 
        
        # 3. Calculer l'angle entre les deux vecteurs
        dot_product = np.dot(rotation_axis, line_of_sight)
        # S'assurer que le produit scalaire est dans [-1, 1]
        dot_product = np.clip(dot_product, -1.0, 1.0)
        aspect_angle_rad = np.arccos(dot_product)
        aspect_angle_deg = np.degrees(aspect_angle_rad)
        
        
        return aspect_angle_deg
    
    def calculate_silhouette_center(self, projected_vertices):
        center = np.mean(projected_vertices, axis=0)
        return center

    def rescale_silhouette(self, projected_vertices, scale_factor):
        
        center = self.calculate_silhouette_center(projected_vertices)
        centered_vertices = projected_vertices - center
        scaled_centered_vertices = centered_vertices * scale_factor
        scaled_vertices = scaled_centered_vertices + center
        
        return scaled_vertices
    
    def determine_silhouette_edges_closed_mesh(self, vertices, triangles, face_visibility):
        """
        Détermine les arêtes de silhouette pour un maillage fermé.
    
        Parameters:
        - vertices: numpy.ndarray, shape (n, 3)
            Coordonnées des sommets du modèle.
        - triangles: numpy.ndarray, shape (m, 3)
            Indices des sommets qui définissent chaque triangle.
        - face_visibility: numpy.ndarray, shape (m,)
            Booléen indiquant si la face est frontale (True) ou arrière (False).
    
        Returns:
        - silhouette_edges: list of tuples
            Liste des arêtes qui composent la silhouette.
        """
        # Création d'un dictionnaire pour stocker les faces adjacentes à chaque arête
        edge_face_dict = {}
    
        for i, tri in enumerate(triangles):
            for j in range(3):
                edge = tuple(sorted((tri[j], tri[(j + 1) % 3])))
                if edge in edge_face_dict:
                    edge_face_dict[edge].append(i)
                else:
                    edge_face_dict[edge] = [i]
    
        silhouette_edges = []
    
        for edge, faces in edge_face_dict.items():
            if len(faces) == 2:
                f1, f2 = faces
                vis1 = face_visibility[f1]
                vis2 = face_visibility[f2]
                if vis1 != vis2:
                    silhouette_edges.append(edge)
            else:
                # Si l'arête n'est pas partagée par deux faces, elle est automatiquement sur la silhouette
                silhouette_edges.append(edge)
    
        return silhouette_edges
    
    
    def calculate_face_normals_and_visibility(self, vertices, triangles, view_direction):
        """
        Calcule les normales des faces et détermine si elles sont frontales ou arrière.
    
        Parameters:
        - vertices: numpy.ndarray, shape (n, 3)
            Coordonnées des sommets du modèle (dans l'espace 3D après projection).
        - triangles: numpy.ndarray, shape (m, 3)
            Indices des sommets qui définissent chaque triangle.
        - view_direction: numpy.ndarray, shape (3,)
            Direction de vue (vecteur).
    
        Returns:
        - face_normals: numpy.ndarray, shape (m, 3)
            Normales des faces.
        - face_visibility: numpy.ndarray, shape (m,)
            Booléen indiquant si la face est frontale (True) ou arrière (False).
        """
        # Calcul des vecteurs des arêtes du triangle
        v0 = vertices[triangles[:, 0]]
        v1 = vertices[triangles[:, 1]]
        v2 = vertices[triangles[:, 2]]
    
        # Vecteurs des côtés du triangle
        edge1 = v1 - v0
        edge2 = v2 - v0
    
        # Calcul des normales des faces
        face_normals = np.cross(edge2, edge1)
    
        # Normalisation des normales
        face_normals /= np.linalg.norm(face_normals, axis=1)[:, np.newaxis]
    
        # Produit scalaire entre la normale et la direction de vue
        dot_products = np.dot(face_normals, view_direction)
    
        # Détermination de la visibilité des faces
        face_visibility = dot_products > 0  # Faces frontales si l'angle est obtus
    
        return face_normals, face_visibility
    
    def __rot_mat(self, alpha, delta):
       
        delta = -np.pi/2 + delta
        alpha *= -1
        
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosd = np.cos(delta)
        sind = np.sin(delta)
        
        R_mat = np.array([[ cosd * cosa, -cosd * sina, sind ],
                          [ sina, cosa, 0 ],
                          [-sind * cosa, sind * sina, cosd]])
        
        
        return R_mat
            
    def project_s(self, model, alpha, delta):
    
        R = self.__rot_mat(alpha, delta)
        
        #vertices_transformed
        if self.shouldPorject:
            vertices = model.vertices_transformed @ R
        else:
            vertices = model.vertices @ R 
        
        return vertices
    
    def project_vertices_on_plane(self, vertices, occ):
        x_eps = np.dot(vertices, occ.sEps)
        y_nu = np.dot(vertices, occ.sNu)
        projected_points = np.column_stack((x_eps, y_nu))
        return projected_points
        
    def compute_projection(self, lamb, beta, angle, model, occ):
        
        view_direction = self.view_dir(occ)
        model.rotate_along_spin(angle + lamb, self.shouldPorject)
        alpha, delta = self.__compute_geocentric_projection(lamb, beta, occ.event.date)
        vertices = self.project_s(model, alpha, delta)
        projected_vertices = self.project_vertices_on_plane(vertices,  occ)
        face_normals, face_visibility = self.calculate_face_normals_and_visibility( vertices, model.triangles, view_direction )
        edge = self.determine_silhouette_edges_closed_mesh( vertices, model.triangles, face_visibility )
        #projected_vertices = self.project_vertices_on_plane(vertices, occ)
    
        return vertices, projected_vertices, edge 
    
    def ecplitic2Equatorial(self, vert, obliquity_time):
        epsilon = np.deg2rad( self.__precise_obliquity(obliquity_time) )

        # Définir la matrice de rotation autour de x
        R_x = np.array([
            [1, 0, 0],
            [0, np.cos(epsilon), -np.sin(epsilon)],
            [0, np.sin(epsilon), np.cos(epsilon)]
        ])
        
        # Appliquer la rotation à tous les vertices
        # On suppose que vertices est un tableau de forme (N, 3)
        vert = vert @ R_x.T  # On transpose car les vertices sont en lignes
        return vert
    
    def compute_projection2(self, lamb, beta, t, period, model, occ):
        
        view_direction = self.view_dir(occ)
        
        vertices = model.in_ecliptic_cord(lamb, beta, t, period, self.shouldPorject)
        vertices = self.ecplitic2Equatorial(vertices, occ.event.date)
        projected_vertices = self.project_vertices_on_plane(vertices,  occ)
        face_normals, face_visibility = self.calculate_face_normals_and_visibility( vertices, model.triangles, view_direction )
        edge = self.determine_silhouette_edges_closed_mesh( vertices, model.triangles, face_visibility )
    
        return vertices, projected_vertices, edge 
    
    def __precise_obliquity(self, date):
        
        t = Time(date, format = "iso")
  

        julian_centuries = (t.jd - 2451545.0) / 36525.0
        obliquity_deg = 23.439292 - 0.0130042 * julian_centuries

        return obliquity_deg
        
    
    def __compute_geocentric_projection(self, lamb, beta, obliquity_time):
       
       lambda_ecl = np.deg2rad(lamb)
       beta_ecl = np.deg2rad(beta)
       
     
       epsilon = np.deg2rad(self.__precise_obliquity(obliquity_time) )# Obliquité de l'écliptique
 

       alpha = np.arctan2(
           np.sin(lambda_ecl) * np.cos(epsilon) - np.tan(beta_ecl) * np.sin(epsilon),
           np.cos(lambda_ecl)
       )
       delta = np.arcsin(
           np.sin(beta_ecl) * np.cos(epsilon) + np.cos(beta_ecl) * np.sin(epsilon) * np.sin(lambda_ecl)
       )

       return alpha, delta
    
    def __compute_avg_size_and_centroide(self, occ, size_factor = 1):
        size = self.__compute_size(occ)
        centroide = self.__compute_centroide(occ)
        
        
        return size, centroide 
             
             
    def __compute_size_and_centroide(self, occ, isSolutionOne =True):
       if isNone(self.initialCentroids) or occ.date().split(" ")[0] not in self.initialCentroids:
            c = self.__compute_centroide(occ)
       else:
            c = np.asarray(self.initialCentroids[occ.date().split(" ")[0]])[::-1]            
       
       solIdx = 0 if isSolutionOne else 1

       if isNone(self.initialSizes) or occ.date().split(" ")[0] not in self.initialSizes or isNone(self.initialSizes[occ.date().split(" ")[0]][solIdx]):
            s = self.__compute_size(occ, c)
       else:
            s = self.initialSizes[occ.date().split(" ")[0]][solIdx] 
      
       return s, c
            
    def __compute_centroide(self, occ):
        
       avgPos = np.zeros((2))
       n = 0

       maxSize = np.zeros((2))

       for cordId in occ.cords:
           cord = occ.cords[cordId]
           status = occ.cordsStatue[cordId]
            
           if status == 0:
               dst = cord[1, :] - cord[0, :]
               if np.linalg.norm(maxSize) < np.linalg.norm(dst):
                   maxSize = dst
               
       for cordId in occ.cords:
           cord = occ.cords[cordId]
           status = occ.cordsStatue[cordId]
            
           if status == 0:
               dst = cord[1, :] - cord[0, :]
               n += 1
               avgPos += (cord[1, :] + cord[0, :] ) / 2
        
       return avgPos / n
                
    def __compute_size(self, occ, centroide = None):

       avgSize = 0
       n = 0
       
       if isNone(centroide):
           centroide = self.__compute_centroide(occ)
       
       for cordId in occ.cords:
           cord = occ.cords[cordId]
           status = occ.cordsStatue[cordId]
           
           if status == 0:
               n += 1
               
               avgSize += np.linalg.norm(centroide - cord[1, : ])
               avgSize += np.linalg.norm(centroide - cord[0, : ])
            
       return avgSize / (2*n)
    
    
    def __setAngle(self, t, model):
        model.solution["angle"] = model.angOfRot(t)
    
    def __optimise(self, initial_param, model, occ, occIdx, nChord, optimiseFlag = 0):
        
        t = occ.t_ref()
    
        #print("point =============", initial_param)
        shouldReproject = False
        
        if model.solution["shouldLambdaBeOptimised"]:
            shouldReproject = True
            lamb = initial_param[0]
        else:
            lamb = model.solution["lambda"]
           
        
        if model.solution["shouldBetaBeOptimised"]:
            shouldReproject = True
            beta = initial_param[1]
        else:
            beta = model.solution["beta"]
            
            
        if model.solution["shouldPeriodBeOptimised"]:
            model.solution["period"] = initial_param[2]
            shouldReproject = True
        
        size = initial_param[3]
        centx = initial_param[4]
        centy = initial_param[5]

        self.__setAngle(t, model)
        #model.rotate_along_spin(model.solution["angle"] + lamb, self.shouldPorject)
        
        #vertices, projected_vertices, edges = self.compute_projection(lamb, beta, model.angOfRot(t), model, occ)
        

        if shouldReproject:
            self.verticesMem[occIdx], self.projectedVertMem[occIdx], self.edgesMem[occIdx] = self.compute_projection2(lamb, beta, t, model.solution["period"], model, occ)

        projected_vertices = self.rescale_silhouette(self.projectedVertMem[occIdx], size)
        projected_vertices += np.asarray([centx, centy])

    
        x2 = self.compute_chi_square( projected_vertices, self.edgesMem[occIdx], occ, nChord, optimiseFlag, model.id)
        
        
        return x2
    
    def __optimiseDependently(self, initial_param, models, nChord, optimiseFlag = 0):
        
        x2 = 0

        for i, occ in enumerate(self.occults):
            
             param = [ initial_param[0], initial_param[1], initial_param[2], initial_param[3], initial_param[ 4 + i * 2 ], initial_param[ 4 + i*2 + 1 ] ]
             
             x = self.__optimise(param, models[i], occ, i, len(occ.cords), optimiseFlag) 
             
             x2 += x / len(self.occults)  
             
             
             # if models[0].id == self.model2.id and i == len(self.occults)-1:
             #     print("step", initial_param, x2, "\n")
             
        return x2
    
    
    def __minimize(self, func, init, bound, args, shouldOptimise, steps = None, axisOfMovement = {}, n = 2, cpt_limit = 50):
        #axisOfMovement est un dic ou les clé correspond a l'indice de l'occ et les valeur doivent être un array 4x4 ou la première ligne est le vecteur normalisé du premier axe et la deuxième ligne du deuxième axe 

        def isBelowUpBound(idx, keys, param, bound):
            
            if isNone(bound[idx][1]):
                return True
            
            if not isIdxFirstInNewAxis(idx, keys) and not isIdxSecondInNewAxis( idx, keys ):
                return param[0] <= bound[idx][1]
            elif isIdxFirstInNewAxis(idx, keys):
                return param[0] <= bound[idx][1] and param[1] <= bound[idx + 1][1]
            elif isIdxSecondInNewAxis( idx, keys ):
                return param[0] <= bound[idx-1][1] and param[1] <= bound[idx][1]
        
        def isUpperLowBound(idx, keys, param, bound):
            
            if isNone(bound[idx][0]):
                return True
            
            if not isIdxFirstInNewAxis(idx, keys) and not isIdxSecondInNewAxis( idx, keys ):
                return param[0] >= bound[idx][0]
            elif isIdxFirstInNewAxis(idx, keys):
                return param[0] >= bound[idx][0] and param[1] >= bound[idx + 1][0]
            elif isIdxSecondInNewAxis( idx, keys ):
                return param[0] >= bound[idx-1][0] and param[1] >= bound[idx][0]
            
        def makeStep(params, step, direction):
            return np.asarray(params) + step * np.asarray(direction)
        
        def isIdxFirstInNewAxis(i, keys):
            return i in np.asarray(keys)*2 + 4 
        
        def isIdxSecondInNewAxis(i, keys):
            return i in np.asarray(keys)*2 + 5 
        
        def binParam(idx, keys, params):
            if not isIdxFirstInNewAxis(idx, keys) and not isIdxSecondInNewAxis( idx, keys ):
                param = [params[idx]]
            elif isIdxFirstInNewAxis(idx, keys):
                param = [params[idx], params[idx + 1]]
            elif isIdxSecondInNewAxis( idx, keys ):
                param = [params[idx - 1], params[idx]]
            else:
                print("error in paramBinning for solver.")
                sys.exit()
                
            return param
        
        def idx2key(idx, keys):
            
            if isIdxFirstInNewAxis(idx, keys):
                return (idx - 4) / 2
            if isIdxSecondInNewAxis(idx, keys):
                return (idx - 5) / 2
            
            return idx
        
        def unBinParam(idx, params, param, keys):
             if not isIdxFirstInNewAxis(idx, keys) and not isIdxSecondInNewAxis(idx, keys ):
                params[idx] = param[0]
             elif isIdxFirstInNewAxis(idx, keys):
                 params[idx] = param[0]
                 params[idx + 1] = param[1]
             elif isIdxSecondInNewAxis( idx, keys ):
                 params[idx - 1] = param[0]
                 params[idx] = param[1]
                 
             else:
                 print('problem ubinning params')
                 sys.exit()
                 
        def makeAMove(i, keys, params, param, step, directions):
            
            if len(param) == 1:
                stepParam = makeStep(param, step, directions)
            else:
                idxDir = 0 if isIdxFirstInNewAxis(i, keys) else 1
                stepParam = makeStep(param, step, directions[idxDir])
             
            unBinParam(i, params, stepParam, keys)

        
        def opti(func, args, shouldOptimise, x2, params, step, bound, axisOfMovement):
            
            keys = list(axisOfMovement.keys())
            
            for i in range(len(params)):
                
                if not shouldOptimise[i]:
                    continue
                
                param = binParam(i, keys, params)

                if len(param) == 1:
                    direction = [1]
                else:
                    direction = axisOfMovement[ idx2key(i, keys) ]
                
                makeAMove( i, keys, params, param, step[i], direction )
                                
                
                if isBelowUpBound(i, keys, param, bound):
                    x2p = func(params, *args)
                else:
                    x2p = x2 + 1
   
                makeAMove(i, keys, params, param, -step[i], direction)
            
                if isUpperLowBound(i, keys, param, bound):
                    x2n = func(params, *args)
                else:
                    x2n = x2 + 1

                if x2p < x2 and x2p <= x2n:
                    makeAMove(i, keys, params, param, step[i], direction)
                    x2 = x2p
                elif x2n < x2 and x2n < x2p:
                    makeAMove(i, keys, params, param, -step[i], direction)
                    x2 = x2n
                else:
                    unBinParam(i, params, param, keys)
                
                
                
            return x2
    
    
        def oneRun(func, args, shouldOptimise, x2i, params, step, bounds, limit, axisOfMovement):

            cpt_good_solution = 0
            cpt = 0
            
            while cpt_good_solution < 4:
                            
                if cpt_good_solution >= 2:
                    x2tpm = opti(func, args, shouldOptimise, x2i, params, step, bound, axisOfMovement)
                else:
                    x2tpm = opti(func, args, shouldOptimise, x2i, params, step, bound, {}) 
                
                if x2tpm == x2i:
                    cpt_good_solution += 1
                else:
                    cpt_good_solution = 0
                
                x2i = x2tpm
                        
                if cpt == limit:
                    print("limit reached. this sub division stop ! ")
                    return x2i
                    
                cpt +=1

            return x2i
        
        if isNone(steps):
            steps = [1, 1, 0.00001, 1]
            
            for i in range(len(init) - 4):
                steps += [1]
        

        init = np.asarray(init)
        steps = np.asarray(steps)
      
        x2 = func(init, *args)
        
        for i in range(n):
            print("subdivisions progress: ", (i / n)*100 , "%")
            x2 = oneRun(func, args, shouldOptimise, x2, init, steps / (i+1), bound, cpt_limit, axisOfMovement)
        
        return init, x2
            

    def compute_chi_square(self, projected_vertices, silhouette_edges, occ, nChord, optimiseFlag = 0, modelid = -1): #optimiseFlag indique si l'optimisation doit se faire sur la valeur de la corde, ou sur les incertitudes: 0 - cord, 1 - incertidudes inférieur, 2 - incertidudes supérieur 
        
    
        def pair_points(P1, P2, Q1, Q2):
            """
            Associe (P_gauche, Q_gauche) et (P_droite, Q_droite)
            en se basant *uniquement* sur la comparaison des x :
              - "plus petit x" = plus à gauche,
              - "plus grand x" = plus à droite.
            En cas d'égalité stricte sur x, on compare y pour trancher.
            
            Retourne la liste :
              [
                 (P_le_plus_a_gauche, Q_le_plus_a_gauche),
                 (P_le_plus_a_droite, Q_le_plus_a_droite)
              ]
            """
            # On convertit en np.array pour être sûr de pouvoir indexer [0], [1]
            P1 = np.array(P1, dtype=float)
            P2 = np.array(P2, dtype=float)
            Q1 = np.array(Q1, dtype=float)
            Q2 = np.array(Q2, dtype=float)
           
            # Petite fonction pour déterminer (ptGauche, ptDroite)
            # selon x (puis y en cas d'ex-æquo).
            def left_right(ptA, ptB):
                if ptA[0] < ptB[0]:
                    return ptA, ptB
                elif ptA[0] > ptB[0]:
                    return ptB, ptA
                else:
                    # x identiques, on tranche sur y
                    if ptA[1] <= ptB[1]:
                        return ptA, ptB
                    else:
                        return ptB, ptA
        
            # On trie P1, P2 "de gauche à droite"
            P_left, P_right = left_right(P1, P2)
        
            # On trie Q1, Q2 "de gauche à droite"
            Q_left, Q_right = left_right(Q1, Q2)
        
            # On associe "gauche avec gauche" et "droite avec droite"
            return [
                (P_left,  Q_left),
                (P_right, Q_right)
            ]
        
        x2 = 0
        
        polygone = Polygon(projected_vertices, silhouette_edges)
        
        meanv = np.mean(projected_vertices, axis = 0)
        maxDist = np.max(np.linalg.norm(projected_vertices[0] - projected_vertices, axis = 1))
        
        for cordId in occ.cords:
            
            status = occ.cordsStatue[cordId]
            weight = occ.cordWeight[cordId]
            err = occ.cordsErr[cordId]
            
            if optimiseFlag == 0 or optimiseFlag == 3 or status == 1:
                cord = occ.cords[cordId]
            
            elif optimiseFlag == 1 and status != 1:
                cord = np.asarray([err["d_err_p"], err["r_err_n"] ])
                    
            elif optimiseFlag == 2 and status != 1:
                cord = np.asarray([err["d_err_n"], err["r_err_p"] ])
    
            if status == 2:
                
                if Event.isDmissing(occ.event[cordId]):
                    point = cord[1]
                    point_far = point + occ.mean_direction * -1e1 
                    line = np.asarray([point_far, point])
                    
                if Event.isRmissing(occ.event[cordId]):
                    point = cord[0]
                    point_far = point + occ.mean_direction * 1e1
                    line = np.asarray([point, point_far])


                intersections = polygone.intersect(line[0], line[1])
                
                if len(intersections) == 2:
                    pair1, pair2 = pair_points(intersections[0], intersections[1], np.array(point)[::-1], np.array(point_far)[::-1])
                    
                    if Event.isDmissing( occ.event[cordId] ):
                       
                        if optimiseFlag == 0:
                            sigma = np.linalg.norm(err["r_err_n"] - err["r_err_p"])
                        else:
                            sigma = 1
                            
                        if sigma == 0:
                            sigma = 1
                        
                        x2 += ( np.linalg.norm(pair2[0] - pair2[1])**2 / sigma**2   ) * weight
                        
                        
                  
                    if Event.isRmissing( occ.event[cordId] ):
                        
                        if optimiseFlag == 0:
                            sigma = np.linalg.norm(err["d_err_n"] - err["d_err_p"])
                        else:
                            sigma = 1
                            
                        if sigma == 0:
                            sigma = 1
                        
                        x2 += ( np.linalg.norm(pair1[0] - pair1[1])**2 / sigma**2 ) * weight
                             
                        
                else:
                   x2 += (1e8 + (np.linalg.norm(meanv - point)**2 - maxDist**2) * nChord ) * weight
                    
            elif status == 0:
                
                intersections = polygone.intersect(cord[0], cord[1]) 
                
                if len(intersections) == 0:
                    
                    meanc = np.mean(cord, axis = 0)
                    
                    x2 += ( 1e8 + (np.linalg.norm(meanv - meanc)**2  - maxDist**2)* nChord ) * weight
             
            
                elif len(intersections) == 2:
                    
                    
                    if optimiseFlag == 0:
                        sigma1 = np.linalg.norm(err["d_err_p"] - err["d_err_n"]) 
                        sigma2 = np.linalg.norm(err["r_err_p"] - err["r_err_n"])
                    else:
                        sigma1 = 1
                        sigma2 = 1
  
        
                    if sigma1 == 0:
                        sigma1 = 1
                    if sigma2 == 0:
                        sigma2 = 1
                    
                    pair1, pair2 = pair_points(intersections[0], intersections[1], np.array(cord[0])[::-1], np.array(cord[1])[::-1])

                    d1 = np.linalg.norm(pair1[0] - pair1[1]) 
                    d2 = np.linalg.norm(pair2[0] - pair2[1]) 
                    
                    d1b = np.linalg.norm(pair1[0] - pair2[1])
                    d2b = np.linalg.norm(pair2[0] - pair1[1])
                    
                    x2 += ( (d1*d1 / sigma1**2 + d2*d2 / sigma2**2) ) * weight
                    #x2 += (d1**2/sigma1**2 + d2**2/sigma2**2 - (d1/sigma1**2 + d2/sigma2**2)**2/(1/sigma1**2 + 1/sigma2**2)) * weight
                    # if modelid == self.model2.id and occ.date().split(" ")[0] == "2024-09-06":
                    #     print("cordid:" , cordId, ( (d1*d1 / sigma1**2 + d2*d2 / sigma2**2) ) * weight, ( (d1*d1 / sigma1**2 + d2*d2 / sigma2**2) / (1 / sigma1**2 + 1 / sigma2**2) ) * weight, pair1, pair2)
                    
   


            elif status == 1:
                # Observation négative: la corde ne doit pas intersecter la silhouette
                # Construire une ligne avec la direction moyenne passant par le point
                point = cord[0]  # cord[0] == cord[1]
                point_far = point + occ.mean_direction * 1e1  # Point éloigné dans la direction moyenne
                negative_line = np.asarray([point - occ.mean_direction * 1e5, point_far])
    
                # Calculer l'intersection avec la silhouette
                intersection = polygone.intersect(negative_line[0], negative_line[1])
                
                if len(intersection) != 0:
                    # Pénaliser si la ligne intersecte la silhouette
                    # La pénalité peut être proportionnelle à la longueur de l'intersection
                    
                    meanInt = np.mean(np.asarray(intersection), axis = 0)
                    dist = np.linalg.norm(meanInt - meanv[::-1])
                    #print(occ.date(), maxDist, dist, maxDist**2  - dist**2)
                    x2 += (1e8  + (maxDist**2  - dist**2) * (1+nChord)) * weight

        return x2
    
    def view_dir(self, occ):
        view_direction = -np.cross(occ.sEps, occ.sNu)
        return view_direction / np.linalg.norm(view_direction)   
    
    def feedParamToSol(self, dicc, solution, isIndependante = True):
        dicc["lambda"] = solution[0]
        dicc["beta"] = solution[1]
        dicc["period"] = solution[2]
        dicc["size"] = solution[3]
        dicc["centroide"] = np.zeros((2))
        dicc["centroide"][0] = solution[4]
        dicc["centroide"][1] = solution[5]
        
    def getParamFromSol(self, dicc):
        
        l = dicc["lambda"]
        b = dicc["beta"]
        p = dicc["period"]
        s = None
        cx = None
        cy = None
        
        if "size" in dicc:
            s = dicc["size"]
        if "centroide" in dicc:
            cx = dicc['centroide'][0]
            cy = dicc["centroide"][1]
        
        return l, b, p, s, cx, cy
            
        
    def __solve(self, func, init, args, shouldOptimise, bounds, solverType, solverStrategy, n_sub_division, sub_max_stop):
        
        def defSteps(init):
            steps = [0.5, 0.5, 0.00001, 0.2]
            
            for i in range(len(init) - 4):
                steps += [0.5]
            
            return np.asarray(steps)
        
        def initListOfParams(param, bounds):
            multiInit = [param.copy() for _ in range((len(param) - 1) * 2 + 1)]
            
            l = param[0].copy()
            b = param[1].copy()
            s = param[3].copy()
            
            for i in range(len(multiInit)):
                multiInit[i][3] = s * 1.10
            
       
            multiInit[0][0] = l + np.abs(bounds[0][1] - l) / 2
            multiInit[1][0] = l - np.abs(bounds[0][0] - l) / 2
 
            multiInit[2][1] = b + np.abs(bounds[1][1] - b) / 2
            multiInit[3][1] = b - np.abs(bounds[1][0] - b) / 2
            multiInit[5][3] = s * 0.90
     
            
            for i in range(int((len(param) - 4) /2)):
                multiInit[6 + i][4 + 2 * i] = param[4 + 2 * i] + 10  # Prend 4, 6, 8, ...
                multiInit[7 + i][5 + 2 * i] = param[5 + 2 * i] - 10  # Prend 5, 7, 9, ...
            
            return multiInit
        
        axisOfPreference = {}
        
        for i, occ in enumerate(self.occults):
            axisOfPreference[i] = np.asarray([self.normVecFromCoef(occ.a), self.normVecFromCoef(-1/occ.a)])
        
       
        
        def updateMult(multi, index):
            index2 = int(index/2)
            if index2 >1:
                index2 += 1
            
            for i in range(len(multi)):
                multi[i][index2] = multi[index][index2]
            print("update", index, index2, multi[index][index2], multi)
            return multi
        
        if solverStrategy == "scipy":
            result = minimize(func, init, args=args, method=solverType, bounds = bounds)
            x2 = func(result.x, *args)
            return result.x, x2
        
        elif solverStrategy == "home":
            return self.__minimize(func, init, bounds, args, shouldOptimise, steps = None, axisOfMovement = axisOfPreference, n = n_sub_division, cpt_limit = sub_max_stop)
        
        elif solverStrategy == "hh":
            result, x2 = self.__minimize(func, init, bounds, args, shouldOptimise, steps = defSteps(init)*10, axisOfMovement = axisOfPreference,  n = n_sub_division, cpt_limit = sub_max_stop)
            result[3] *= 1.15
            return self.__minimize(func, result, bounds, args, shouldOptimise, steps = None, axisOfMovement = axisOfPreference, n = n_sub_division, cpt_limit = sub_max_stop)
            
        elif solverStrategy == "hybrid":
            result = minimize(func, init, args=args, method=solverType, bounds = bounds)
            result.x[3] *= 1.1
            return self.__minimize(func, result.x, bounds, args, shouldOptimise, steps = defSteps(init), axisOfMovement = axisOfPreference, n = n_sub_division, cpt_limit = sub_max_stop)
        
        elif solverStrategy == "hybrid_r":
            res, _ = self.__minimize(func, init, bounds, args, shouldOptimise, steps = None, axisOfMovement = axisOfPreference, n = n_sub_division, cpt_limit = sub_max_stop)
            res[3] *= 1.1
            result = minimize(func, res, args=args, method=solverType, bounds = bounds)
            x2 = func(result.x,  *args)
            
            return result.x, x2
        
        elif solverStrategy == "deep_hybrid":
            result = minimize(func, init, args=args, method=solverType, bounds = bounds)
           
            x2 = func(result.x, *args)
            print("initial x2: ", x2)
            finalParam = result.x 
            
            multiInit = initListOfParams(result.x, bounds)
            steps = defSteps(init)

            for i in range(len(multiInit)):
                print("progress:", (i / len(multiInit)) * 100, "%")

                res, x2tpm = self.__minimize(func, multiInit[i], bounds, args, shouldOptimise, steps = steps, axisOfMovement = axisOfPreference, n = n_sub_division, cpt_limit = sub_max_stop)
                
                if x2tpm < x2:
                    x2 = x2tpm
                    finalParam = res
                
            return finalParam, x2
        
        elif solverStrategy == "deep_hybrid_r":
            finalParam, x2 = self.__minimize(func, init, bounds, args, shouldOptimise, steps = None, axisOfMovement = axisOfPreference, n = n_sub_division, cpt_limit = sub_max_stop)
            
            multiInit = initListOfParams(result.x, bounds)
            
            for i, initParam in enumerate(multiInit):
                print("progress:", (i / len(multiInit)) * 100, "%")
                result = minimize(func, initParam, args=args, method=solverType, bounds = bounds)
                x2tpm = func(result.x, *args)
                
                if x2tpm < x2:
                    x2 = x2tpm
                    finalParam = result.x
                    
                    
            return finalParam, x2
        
        else:
            print("The strategy ", solverStrategy," do not exist. Please select one of this: scipy, home, hybrid, deep_hybrid, hybrid_r, deep_hybrid_r")
            sys.exit()
                                        

    def determinMajorMinorAxisOrSilhouette(self, pv):
        
        maxDist = 0
        p1Max = np.zeros((2))
        p2Max = np.zeros((2))
        
        for p1 in pv:
            for p2 in pv:
                dist = np.linalg.norm(p2 - p1)
                if dist > maxDist:
                    maxDist = dist
                    p1Max = p1
                    p2Max = p2
            
        major_vector = p2Max - p1Max
        major_norm = major_vector / np.linalg.norm(major_vector)
        perp_vector = np.array([-major_norm[1], major_norm[0]])
        
       
        pv = np.array(pv)
        c_coords = pv.dot(major_norm)
        d_coords = pv.dot(perp_vector)
        
        
        minor_axis = d_coords.max() - d_coords.min() 
        
        tol = (c_coords.max() - c_coords.min()) / 100.0  
        
        max_width = 0  
        best_pair = None
        n = len(pv)
        for i in range(n):
            for j in range(i+1, n):
                if abs(c_coords[i] - c_coords[j]) < tol:
                    width = abs(d_coords[i] - d_coords[j])
                    if width > max_width:
                        max_width = width
                        best_pair = (i, j)
                        
        minor_axis = max_width
        
        majAxis = {"dist" : maxDist,
                   "p1" : p1Max,
                   "p2" : p2Max}
        
        minAxis = {"dist" : minor_axis,
                   "p1": pv[best_pair[0]],
                   "p2": pv[best_pair[1]]}
        
        return majAxis, minAxis
                                          
    
    def findSolutionOnCircle(self, x, y, a, b, s):
        
        a2nd = (1 + a*a)
        b2nd = -2*(x + y*a - a*b)
        c2nd = x*x - s*s + (y - b)*(y - b)
        
        delt = b2nd*b2nd - 4 * a2nd * c2nd
        
        return (- b2nd - np.sqrt(delt)) / (2*a2nd), (-b2nd + np.sqrt(delt)) / (2*a2nd)
    
    def normVecFromCoef(self, a):
        return np.asarray( [a / np.sqrt(1 + a*a), 1 / np.sqrt(1 + a*a)] )
    
    def multiInitPoints(self, initialPoint, model, occ, idxOfOcc = 0):
        
        lamb = initialPoint[0]
        beta = initialPoint[1]
        period = initialPoint[2]
        size = initialPoint[3]
        
        cy = initialPoint[ 4 + idxOfOcc * 2 ] 
        cx = initialPoint[ 5 + idxOfOcc * 2 ]
        
        v, pv, e = self.compute_projection2(lamb, beta, occ.t_ref(), period, model, occ)
        
        major, minor = self.determinMajorMinorAxisOrSilhouette(pv)
        
        print("major", major, "minor", minor)
        
        init = [initialPoint.copy() for i in range(9)]
        
        amajor = (major["p2"][0] - major["p1"][0]) / (major["p2"][1] - major["p1"][1])
        a = -1/occ.a
        b = cy - a*cx
        
        theta = np.arctan(abs((amajor - a) / (1 + a * amajor)))
        
        px1, px2 = self.findSolutionOnCircle(cx, cy, a, b, size *(1 + major["dist"] / minor["dist"] * np.cos(theta)))
        py1, py2 = a*px1 + b, a*px2 + b
        
        print(px1, py1, px2, py2)
        
        ag = -1/a
        bg1 = py1 - ag*px1 
        bg2 = py2 - ag*px2
        b = cy - ag * cx
        
        cx1, cx2 = self.findSolutionOnCircle(px1, py1, ag, bg1, size)
        cx3, cx4 = self.findSolutionOnCircle(px2, py2, ag, bg2, size)
 
        cy1, cy2, cy3, cy4 = ag * cx1 + bg1, ag * cx2 + bg1, ag * cx3 + bg2, ag * cx4 + bg2
        
        px3, px4 = self.findSolutionOnCircle(cx, cy, ag, b, size *(1 + major["dist"] / minor["dist"] * np.cos(theta)))
        py3, py4 = ag*px3 + b, ag*px4 + b
 
        init[1][4 + idxOfOcc*2 + 1] = cx1 
        init[1][4 + idxOfOcc*2] = cy1 
        
        init[2][4 + idxOfOcc*2 + 1] = cx2 
        init[2][4 + idxOfOcc*2] = cy2
        
        init[3][4 + idxOfOcc*2 + 1] = cx3 
        init[3][4 + idxOfOcc*2] = cy3 
        
        init[4][4 + idxOfOcc*2 + 1] = cx4 
        init[4][4 + idxOfOcc*2] = cy4
        
        init[5][4 + idxOfOcc*2 + 1] = px1 
        init[5][4 + idxOfOcc*2] = py1
        
        init[6][4 + idxOfOcc*2 + 1] = px2 
        init[6][4 + idxOfOcc*2] = py2
        
        init[7][4 + idxOfOcc*2 + 1] = px3 
        init[7][4 + idxOfOcc*2] = py3
        
        init[8][4 + idxOfOcc*2 + 1] = px4 
        init[8][4 + idxOfOcc*2] = py4
        
        
        return init
        
    
    def occShouldntBeMultiInit(self, occ, isSolution1 = True):
        
        sol = 1 if isSolution1 else 2
        
        return occ.date().split(" ")[0] not in self.multi or (sol not in self.multi[occ.date().split(" ")[0]] and 3 not in self.multi[occ.date().split(" ")[0]]) # 1 - solution1, 2 - solution2, 3 - solution 1 and 2
        
                                
    def multiInitialisation(self, initial_param, model, occ, occIdx, nChords, shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop, shouldFitDependently = True):
        
        
        
        if shouldFitDependently:
            initial_params = self.multiInitPoints(initial_param, model[occIdx], occ, occIdx)
            i = occIdx
        else:
            initial_params = self.multiInitPoints(initial_param, model, occ, 0)
            i = 0
                
        x2 = [0] * len(initial_params)
        results = [None] * len(initial_params)
            
        for j, init in enumerate(initial_params):
           
            if shouldFitDependently:
                results[j], x2[j] = self.__solve(self.__optimiseDependently, init, (model, nChords, 0), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop) 
            else:
                results[j], x2[j] = self.__solve(self.__optimise, init, (model, occ, occIdx, nChords, 0), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop) 
        
        if np.all(np.isnan(x2)):
            print("!!!!!!!!!! WARNING, only nan detected. multi initialisation aborded. If the issue remind, try another optimizer strategy like hh !!!!!!!!!!!!!!")
            return initial_param
        
        idxMin = np.where(np.asarray(np.abs(x2)) == np.nanmin(np.asarray(np.abs(x2))))[0][0]
        best = results[ idxMin ]
             
        initial_param[4 + i*2] = best[4 + i*2]
        initial_param[5 + i*2] = best[5 + i*2]
        
        return initial_param
            
    def initShouldOptimise(self, init_param, model):
        shouldOptimise = [True] * len(init_param)
        shouldOptimise[0] = bool(model.solution["shouldLambdaBeOptimised"])
        shouldOptimise[1] = bool(model.solution["shouldBetaBeOptimised"])
        shouldOptimise[2] = bool(model.solution["shouldPeriodBeOptimised"])
        
        return shouldOptimise
            
    def fitDependently(self, solverType, solverStrat, n_sub_division, sub_max_stop):
        
        
        if self.ignorUncer:
            flag = 3
        else:
            flag = 0
        
        print(flag) 
        
        solverStratErr = solverStrat
                
        if solverStratErr == "deep_hybrid_r":
            solverStratErr = "hybrid_r"
            print('deep_hybrid_r not allow for uncertainties. Automaticly changed by hybrid_r')
        elif solverStratErr == "deep_hybrid":
            solverStratErr = "hybrid"
            print('deep_hybrid not allow for uncertainties. Automaticly changed by hybrid')
        
        
        models1 = [copy.deepcopy(self.model1) for _ in range(len(self.occults))]
        self.model1.solution["size"] = 0
        
        if not isNone(self.model2):
            models2 = [copy.deepcopy(self.model2) for _ in range(len(self.occults))]
            self.model2.solution["size"] = 0  
        
        nChords = 0

        sizeAvg1 = 0
        sizeAvg2 = 0
        
        self.__initMem(self.model1)

        for i, occ in enumerate(self.occults):
             
             nChords += len(occ.cords)
             s, c = self.__compute_size_and_centroide(occ)

             if s > sizeAvg1:
                 sizeAvg1 = s 
                 models1[i].solution["size"] = s
            
             models1[i].solution["centroide"] = c
            
             
             if not isNone(self.model2): 
                 
                 s2, c2 = self.__compute_size_and_centroide(occ, False)
                 
                 if s2 > sizeAvg2:
                     sizeAvg2 = s2 
                     models2[i].solution["size"] = s2
                     
                 
                 models2[i].solution["centroide"] = c2
             
        l1, b1, p1, _, _, _ = self.getParamFromSol( self.model1.solution )
        
        initial_param = [l1, b1, p1, sizeAvg1*1.15]
        bond = [(l1 - self.model1.solution["err_l_m"], l1 + self.model1.solution["err_l_p"]),
                (b1 - self.model1.solution["err_b_m"], b1 + self.model1.solution["err_b_p"]),
                (None, None),
                (0, None)]
        
        for m in models1:
            initial_param += [m.solution["centroide"][0], m.solution["centroide"][1]]
            bond += [(None, None), (None , None)]
        
        shouldOptimise = self.initShouldOptimise(initial_param, self.model1)
 
        for i, occ in enumerate(self.occults):
            
            if self.occShouldntBeMultiInit(occ, True):
                continue
            

            initial_param = self.multiInitialisation(initial_param, models1, occ, i, nChords, shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)
                
                
                
        
        print("=========Optimisation solution 1======")
        print("initial param", initial_param)
        print("optimise solution")
        
        result, x21 = self.__solve(self.__optimiseDependently, initial_param, (models1, nChords, flag), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop) 
        print("optimisation end with x2:", x21)
        print("evaluate uncertainties")
      
        init_uncert = result.copy()
        init_uncert[3] *= 1.05
        
        
        result_n, x21errn = self.__solve(self.__optimiseDependently, init_uncert, (models1, nChords, 1), shouldOptimise, bond, solverType, solverStratErr, n_sub_division, sub_max_stop) 
        result_p, x21errp = self.__solve(self.__optimiseDependently, init_uncert, (models1, nChords, 2), shouldOptimise, bond, solverType, solverStratErr, n_sub_division, sub_max_stop) 
        print("uncertainties determination end with x2(+):", x21errp, "x2(-)", x21errn)
        
        
        
        for i, occ in enumerate(self.occults):
            self.feedParamToSol(models1[i].solution, list(result[:4]) + list(result[(4 + 2*i) : 4 + 2*(i + 1)]))#add x
     
        
        print("========Solution 1 optimised ! =======")
     
         
        
        if not isNone(self.model2): 
            
            self.__initMem(self.model2)
            
            l2, b2, p2, _, _, _ = self.getParamFromSol( self.model2.solution )
            
            initial_param = [l2, b2, p2, sizeAvg2*1.15]
            bond = [(l2 - self.model2.solution["err_l_m"], l2 + self.model2.solution["err_l_p"]),
                    (b2 - self.model2.solution["err_b_m"], b2 + self.model2.solution["err_b_p"]),
                    (p2 - self.model2.solution["err_p_m"], p2 + self.model2.solution["err_p_p"]),
                    (0, None)]
            
            for i, m in enumerate(models2):
                
                print(self.occults[i].date(), m.solution["centroide"])
                initial_param += [m.solution["centroide"][0], m.solution["centroide"][1]]
                bond += [(None, None), (None , None)]
    
    
            shouldOptimise = self.initShouldOptimise(initial_param, self.model2)
            
            
            for i, occ in enumerate(self.occults):
                
                if self.occShouldntBeMultiInit(occ, False):
                    continue
                
                initial_param = self.multiInitialisation(initial_param, models2, occ, i, nChords, shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)
    
            print("=========Optimisation solution 2======")
            print("initial param", initial_param)
            print("optimise solution ") 
            result2, x22 = self.__solve(self.__optimiseDependently, initial_param, (models2, nChords, flag), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop) 
            print("optimisation end with x2:", x22)
            print("evaluate uncertainties")
            
            init_uncert = result2.copy()
            init_uncert[3] *= 1.05
            
            result2_n, x22errn = self.__solve(self.__optimiseDependently, init_uncert, (models2, nChords, 1), shouldOptimise, bond, solverType, solverStratErr, n_sub_division, sub_max_stop) 
            result2_p, x22errp = self.__solve(self.__optimiseDependently, init_uncert, (models2, nChords, 2), shouldOptimise, bond, solverType, solverStratErr, n_sub_division, sub_max_stop) 
            

            print("uncertainties determination end with x2(+):", x22errp, "x2(-)", x22errn)
            
            for i, occ in enumerate(self.occults):
                self.feedParamToSol(models2[i].solution, list(result2[:4]) + list(result2[(4 + 2*i) : 4 + 2*(i + 1)]))
            
            print("========Solution 2 optimised ! =======")
            
            
            #====================display solutions===================
            
        eqV = models1[0].equivalentVolRadius(models1[i].solution["size"]) * 2
        eqVp = np.abs(eqV - models1[i].equivalentVolRadius(result_p[3]) * 2)
        eqVn = np.abs(eqV - models1[i].equivalentVolRadius(result_n[3]) * 2)
        strSol1 = "======================Results============================\n"
        dTxt = r" $D = {eqV}^{+eqVp}_{-eqVn}$".replace("eqVp", str(round(eqVp, 2))).replace("eqVn", str(round(eqVn, 2))).replace("eqV", str(round(eqV, 2)))
        
        if not isNone(self.model2):
            eqV2 = models2[0].equivalentVolRadius(models2[i].solution["size"]) * 2
            eqVp2 = np.abs(eqV - models2[i].equivalentVolRadius(result2_p[3]) * 2)
            eqVn2 = np.abs(eqV - models2[i].equivalentVolRadius(result2_n[3]) * 2)
            
            strSol2 = "---------------------------------------------------\n"
            
            
            dTxt2 = r" $D = {eqV}^{+eqVp}_{-eqVn}$".replace("eqVp", str(round(eqVp2, 2))).replace("eqVn", str(round(eqVn2, 2))).replace("eqV", str(round(eqV2, 2)))
      
        
        for i, occ in enumerate(self.occults):
            
            strSol1 += "DATE : " + str(occ.date()) + '\n'       
            strSol1 += "Solution: " + str(models1[i].solution) + "\n"
            strSol1 += "Aspect angle: " + str(self.compute_aspect_angle( models1[i], occ)) + "\n\n"
            
            if not isNone(self.model2):
                self.plot_fit_point(occ, models1[i].solution, occ.getAst() + ' ' + occ.date().split(" ")[0] + '\n' + r"Shape 1 ${\chi}^{2}$: " + str(round(x21, 2)) + dTxt  , self.model1, models2[i].solution, r"Shape 2 ${\chi}^{2}$: " + str(round(x22, 2)) + dTxt2, self.model2, "final best", block = True)
                strSol2 += "DATE : " + str(occ.date()) + '\n'       
                strSol2 += "Solution: " + str(models2[i].solution) + "\n"
                strSol2 += "Aspect angle: " + str(self.compute_aspect_angle( models2[i], occ)) + "\n\n"
            else:
                self.plot_fit_point(occ, models1[i].solution, occ.getAst() + ' ' + occ.date().split(" ")[0] + '\n' + r"Shape 1 ${\chi}^{2}$: " + str(round(x21, 2)) + dTxt  , self.model1, None, "", self.model2, "final best", block = True)
        
        strSol1 += "equivalent Volume Diameter = " + str(eqV) + " + " + str(eqVp) + " - " + str(eqVn) + '\n' + "chi2: " + str(x21) + "\n"
        print(strSol1)
        
        if not isNone(self.model2):
            strSol2 += "equivalent Volume Diameter = " + str(eqV2) + " + " + str(eqVp2) + " - " + str(eqVn2) + '\n' + "chi2: " + str(x22) + "\n\n"
            print(strSol2)
       
       

    def fitIndependently(self, solverType, solverStrat, n_sub_division, sub_max_stop):
       
       
       if self.ignorUncer:
           flag = 3
       else:
           flag = 0
       

       
       solverStratErr = solverStrat
               
       if solverStratErr == "deep_hybrid_r":
           solverStratErr = "hybrid_r"
           print('deep_hybrid_r not allow for uncertainties. Automaticly changed by hybrid_r')
       elif solverStratErr == "deep_hybrid":
           solverStratErr = "hybrid"
           print('deep_hybrid not allow for uncertainties. Automaticly changed by hybrid') 
       
       
        
       l1, b1, p1, _, _, _ = self.getParamFromSol( self.model1.solution )
       if not isNone(self.model2):
           l2, b2, p2, _, _, _ = self.getParamFromSol( self.model2.solution )
        
       self.__initMem(self.model1)    
       self.__initMem(self.model2)

       for i, occ in enumerate(self.occults):
     
            s1, c1 = self.__compute_size_and_centroide(occ, True)
            
            self.model1.solution["size"] = s1
            self.model1.solution["centroide"] = c1
            
            initial_param = [l1, b1, p1, s1*2, c1[0], c1[1] ]
           
            bond = [(l1 - self.model1.solution["err_l_m"], l1 + self.model1.solution["err_l_p"]),
                    (b1 - self.model1.solution["err_b_m"], b1 + self.model1.solution["err_b_p"]),
                    (p1 - self.model1.solution["err_p_m"], p1 + self.model1.solution["err_p_p"]),
                    (0, None),
                    (None, None), 
                    (None , None)]
            
            
            shouldOptimise = self.initShouldOptimise(initial_param, self.model1)
            
            if not self.occShouldntBeMultiInit(occ, False):
                initial_param = self.multiInitialisation(initial_param, self.model1, occ, i, len(occ.cords), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop, False)
            
            print("=========Optimisation solution 1======")
            print("initial param", initial_param)
            print("optimise solution")
            result, x21 = self.__solve(self.__optimise, initial_param, (self.model1, occ, i, len(occ.cords), flag), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)
            print("optimisation end with x2:", x21)
            print("evaluate uncertainties")
            init_uncert = result.copy()
            init_uncert[3] *= 1.05
            result_errn, x21errn = self.__solve(self.__optimise, init_uncert, (self.model1, occ, i, len(occ.cords), 1), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)
            result_errp, x21errp = self.__solve(self.__optimise, init_uncert, (self.model1, occ, i, len(occ.cords), 2), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)

            print("uncertainties determination end with x2(+):", x21errp, "x2(-)", x21errn)
        
            if not isNone(self.model2): 
               
                s2, c2 = self.__compute_size_and_centroide(occ, False)
                
                self.model2.solution["size"] = s2
                self.model2.solution["centroide"] = c2
                
                initial_param2 = [l2, b2, p2, s2*2, c2[0], c2[1]]
                
                bond = [(l2 - self.model2.solution["err_l_m"], l2 + self.model2.solution["err_l_p"]),
                        (b2 - self.model2.solution["err_b_m"], b2 + self.model2.solution["err_b_p"]),
                        (p2 - self.model2.solution["err_p_m"], p2 + self.model2.solution["err_p_p"]),
                        (0, None),
                        (None, None), 
                        (None, None)]
                  
                shouldOptimise = self.initShouldOptimise(initial_param, self.model2)  
              
                if not self.occShouldntBeMultiInit(occ, False):
                    initial_param = self.multiInitialisation(initial_param, self.model2, occ, i, len(occ.cords), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop, False)
                
     
                print("=========Optimisation solution 2======")
                print("initial param", initial_param2)
                print("optimise solution")      
                result2, x22 = self.__solve(self.__optimise, initial_param2, (self.model2, occ, i, len(occ.cords), flag), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)
                print("optimisation end with x2:", x22)
                print("evaluate uncertainties")
                init_uncert = result2.copy()
                init_uncert[3] *= 1.05
                result2_errn, x22errn = self.__solve(self.__optimise, init_uncert, (self.model2, occ, i, len(occ.cords), 1), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)
                result2_errp, x22errp = self.__solve(self.__optimise, init_uncert, (self.model2, occ, i, len(occ.cords), 2), shouldOptimise, bond, solverType, solverStrat, n_sub_division, sub_max_stop)
                print("uncertainties determination end with x2(+):", x22errp, "x2(-)", x22errn)

                
                
            
            self.feedParamToSol(self.model1.solution, result)
            
            eqV = self.model1.equivalentVolRadius(self.model1.solution["size"]) * 2
            eqVp = np.abs(eqV - self.model1.equivalentVolRadius(result_errn[3]) * 2)
            eqVn = np.abs(eqV - self.model1.equivalentVolRadius(result_errn[3]) * 2)
        
            dTxt = r" $D = {eqV}^{+eqVp}_{-eqVn}$".replace("eqVp", str(round(eqVp, 2))).replace("eqVn", str(round(eqVn, 2))).replace("eqV", str(round(eqV, 2)))
            
            
        
            if not isNone(self.model2):
                self.feedParamToSol(self.model2.solution, result2)
            
                eqV2 = self.model2.equivalentVolRadius(self.model2.solution["size"]) * 2
                eqVp2 = np.abs(eqV - self.model2.equivalentVolRadius(result2_errp[3]) * 2)
                eqVn2 = np.abs(eqV - self.model2.equivalentVolRadius(result2_errn[3]) * 2)
                
                dTxt2 = r" $D = {eqV}^{+eqVp}_{-eqVn}$".replace("eqVp", str(round(eqVp2, 2))).replace("eqVn", str(round(eqVn2, 2))).replace("eqV", str(round(eqV2, 2)))
            
                self.plot_fit_point(occ, self.model1.solution, occ.getAst() + ' ' + occ.date().split(" ")[0] + '\n' + r"Shape 1 ${\chi}^{2}$: " + str(round(x21, 2)) + dTxt, self.model1, self.model2.solution, r"Shape 2 ${\chi}^{2}$: " + str(round(x22, 2)) + dTxt2, self.model2, "final best", block = True)
            else:
                self.plot_fit_point(occ, self.model1.solution, occ.getAst() + ' ' + occ.date().split(" ")[0] + '\n' + r"Shape 1 ${\chi}^{2}$: " + str(round(x21, 2)) + dTxt, self.model1, None, "", self.model2, "final best", block = True)
            
            print("==================================================")
            print("solution1: ", self.model1.solution)
            print("initial_point:", initial_param)
            print("equivalent diameter = ", eqV, "+", eqVp, "-", eqVn)
            print("Aspect angle = ", self.compute_aspect_angle(self.model1, occ))
            
            if not isNone(self.model2):
                print("-------------------------------------------")
                print("solution2: ", self.model2.solution)
                print("initial_point:", initial_param2)
                print("equivalent diameter = ", eqV2, "+", eqVp2, "-", eqVn2)
                print("Aspect angle = ", self.compute_aspect_angle(self.model2, occ))
                print("\n\n")
            else:
                print("\n\n")
    
    
    def scan(self):
        
        
        def write(param, x2):
            
            nc = int((len(param) - 4) / 2)
            
            strr = str(param[3]) + " "
            
            for i in range(nc):
                strr += str(param[4 + 2*i]) + " " + str(param[5 + 2*i]) + " "
                
            strr += str(x2) + '\n'
            
            return strr
            
            
        
        init = [ 5.39000000e+01, 3.70000000e+00, 2.38963717e+01, 9.18783755e+01, 4.12035511e+03,  4.03976275e+03,  5.58411049e+03, -5.38155180e+01, 2.20495560e+03, -4.21907033e+03]
        
       
        
        cp = 30
        cn = 30
        sp = 1.1
        sn = 0.9
      
        figure_to_adapte = 0
        cx =  np.linspace(init[4 + (figure_to_adapte * 2)] - cn, init[4 + (figure_to_adapte * 2)] + cp, 60)
        cy =  np.linspace(init[5 + (figure_to_adapte * 2)] - cn, init[5 + (figure_to_adapte * 2)] + cp, 60)
        figure_to_adapte = 1
        cx1 = np.linspace(init[4 + (figure_to_adapte * 2)] - cn, init[4 + (figure_to_adapte * 2)] + cp, 60)
        cy1 =  np.linspace(init[5 + (figure_to_adapte * 2)] - cn, init[5 + (figure_to_adapte * 2)] + cp, 60)
        figure_to_adapte = 2
        cx2 = np.linspace(init[4 + (figure_to_adapte * 2)] - cn, init[4 + (figure_to_adapte * 2)] + cp, 60)
        cy2 =  np.linspace(init[5 + (figure_to_adapte * 2)] - cn, init[5 + (figure_to_adapte * 2)] + cp, 60)
       
        size =  np.linspace(init[3] * sn, init[3] * sp, 20)
    
        models = [copy.deepcopy(self.model1) for _ in range(len(self.occults))]
    
        nChords = 0

        for i, occ in enumerate(self.occults):
             
             nChords += len(occ.cords)
             
        f = open("scan.txt", "a")
        for s in size:
        
            figure_to_adapte = 0
           
            for x in cx:
                for y in cy:
                
                    
                    init[4 + (figure_to_adapte * 2)] = x
                    init[5 + (figure_to_adapte * 2)] = y
                    init[3] = s
                    
                    x2 = self.__optimiseDependently(init, models, nChords, 0)
                    
                    strr = write(init, x2)
                    f.write(strr)
          
            figure_to_adapte = 1
            
            for x in cx1:
                for y in cy1:
                
                    
                    init[4 + (figure_to_adapte * 2)] = x
                    init[5 + (figure_to_adapte * 2)] = y
                    init[3] = s
                    
                    x2 = self.__optimiseDependently(init, models, nChords, 0)
                    
                    strr = write(init, x2)
                    
                    f.write(strr)
            
            figure_to_adapte = 1    
            
            for x in cx1:
                for y in cy1:
                 
                     
                     init[4 + (figure_to_adapte * 2)] = x
                     init[5 + (figure_to_adapte * 2)] = y
                     init[3] = s
                     
                     x2 = self.__optimiseDependently(init, models, nChords, 0)
                     
                     strr = write(init, x2)
                     
                     f.write(strr)
                
        f.close()
            

    def fit(self, shouldFitIndependently = False, solverType = "SLSQP", solverStrat = "hybrid", n_sub_division = 2, sub_max_stop = 50, multi = None, initial_centroids = None, intial_size = None, ignorUncertainties = False):
       
        print("strategy selected : ", solverStrat)
        
        self.multi = multi
        self.initialCentroids = initial_centroids
        self.initialSizes = intial_size
        self.ignorUncer = ignorUncertainties
        
        if solverStrat != "home":
            print("solver type selected : ", solverType)
        if solverStrat != "scipy":
            print("number of subdivision:", n_sub_division)
            print("max iteration: ", sub_max_stop)

        if shouldFitIndependently:
           self.fitIndependently(solverType, solverStrat, n_sub_division, sub_max_stop)
        else:
           #self.scan()
           self.fitDependently(solverType, solverStrat, n_sub_division, sub_max_stop)
           

    
    def plot_silhouette_from_vertices_edge(self, projected_vertices, silhouette_edges):
        
        fig, ax = plt.subplots()
        
        for edge in silhouette_edges:
            p1, p2 = edge
            x_values = [projected_vertices[p1][0], projected_vertices[p2][0]]
            y_values = [projected_vertices[p1][1], projected_vertices[p2][1]]
            ax.plot(x_values, y_values, color='blue', linewidth=1)
   
        ax.set_xlabel('n (eta)')
        ax.set_ylabel('ε (xi)')
        ax.set_title('Silhouette de l\'objet projeté sur le plan ε-ν')
        plt.gca().set_aspect('equal', adjustable='box')
        ax.grid(True)
        plt.show(block=False)
    
    def plot_silhouette(self, projected_vertices, edges, c = "blue", ax = None, block = False, centroide = None):
        ax1 = ax
        if ax == None:
            fig, ax1 = plt.subplots(figsize=(8, 8))
     
        label = ""
        
        
        if isNone(centroide):
            centroidex = 0
            centroidey = 0  
            centroide = np.zeros((2))
        else:
            centroidex = centroide[1]
            centroidey = centroide[0]
            centroide = np.asarray(centroide)
            
        for point in projected_vertices:
            point[1] -= centroidex
            point[0] -= centroidey

        
        for i, edge in enumerate(edges):
            p1, p2 = edge
            point1 = projected_vertices[p1]
            point2 = projected_vertices[p2]
            
            if c == "blue" and i == len(edges) - 1:
                label = "shape 1"
            elif i == len(edges) - 1:
                label = "shape 2"
            
            ax1.plot([point1[1], point2[1]], [point1[0], point2[0]], color=c, linewidth=1, label = label)
    
        if ax  == None:
            plt.show(block=block)
    
        return ax1
    
    

    def plot_fit_point(self, occ,  point1, title1, model1, point2 = None, title2 = "", model2 = None, title = "", block = False):
        
        lamb = point1["lambda"]
        beta = point1["beta"]
        ang = point1["angle"]
        centroide = point1["centroide"]
        size = point1["size"]
        
        #vertices1, projected_vertices1, edges1 = self.compute_projection(lamb, beta, ang, model1, occ)
        vertices1, projected_vertices1, edges1 = self.compute_projection2(lamb, beta, occ.t_ref(), model1.solution["period"], model1, occ)
        projected_vertices1 = self.rescale_silhouette(projected_vertices1, size)
        projected_vertices1 += centroide
        
        title = title1
        
        if point2 != None and title2 != "" and model2 != None:
            lamb = point2["lambda"]
            beta = point2["beta"]
            ang = point2["angle"]
            centroide = point2["centroide"]
            size = point2["size"]
            
            #vertices2, projected_vertices2, edges2 = self.compute_projection(lamb, beta, ang, model2, occ)
            vertices2, projected_vertices2, edges2 = self.compute_projection2(lamb, beta, occ.t_ref(), model2.solution["period"], model2, occ)
            projected_vertices2 = self.rescale_silhouette(projected_vertices2, size)
            projected_vertices2 += centroide
        
            title += '\n' + title2
        
        ax = self.plot_fit_cords(occ, projected_vertices1, edges1, title, block = False, centroide = centroide)
       
       
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
       
        if point2 != None and title2 != "" and model2 != None: 
            self.plot_silhouette(projected_vertices2, edges2, c = "m", ax = ax, block = block, centroide=centroide)
          
        ax.set_xlim(xlim[0] - 10, xlim[1] + 10)
        ax.set_ylim(ylim[0] - 10, ylim[1] + 10)
        
        #ax.legend()
        self.plot_projected_spin_axis(ax, model1, model2, occ)        
        ax.title.set_fontsize(16)
        return ax
    
    
    def plot_projected_spin_axis(self, ax, model, model2, occ, position=(0, 0), size=0.08):
        
       
        def addMarker(ax, model, occ, color, label, lw):
            lamb = model.solution['lambda']
            beta = model.solution['beta']
            alpha, delta = self.__compute_geocentric_projection(lamb, beta, occ.event.date)  
            rotation_axis = self.__spherical2cartesian(alpha, delta)  
            
            rotation_axis /=  np.linalg.norm(rotation_axis)
            
           
            proj_x = np.dot(rotation_axis, occ.sNu)
            proj_y = np.dot(rotation_axis, occ.sEps)

            angle = np.radians(self.compute_aspect_angle(model, occ))  
            norm = np.sqrt(proj_x**2 + proj_y**2)
        
            k = 1 / norm
            
            proj_x *= k
            proj_y *= k
            
            marker = 'o' if -np.pi/2 <= angle <= np.pi/2 else 'x'
            s = 10 if marker == "o" else 20
            
            scatter = ax.scatter(proj_x*np.sin(angle), proj_y*np.sin(angle), color=color, marker=marker, s=s, lw=lw, label = label)
            
            #plt.savefig(f"{occ.event.date}.png")
            #plt.savefig(f"{occ.event.date}.eps")
            
            return scatter
        
        inset_ax = ax.inset_axes([position[0], position[1], size, size], transform=ax.transAxes)
        inset_ax.set_xticks([])
        inset_ax.set_yticks([])
        inset_ax.set_frame_on(False)
        circle = plt.Circle((0, 0), 1, color='black', fill=False, lw=1.5)
        inset_ax.add_patch(circle)
        
        existing_legend = ax.get_legend()
        if existing_legend:
            try:
                existing_handles, existing_labels = existing_legend.legendHandles, [text.get_text() for text in existing_legend.get_texts()]
            except:
                existing_handles, existing_labels = existing_legend.legend_handles, [text.get_text() for text in existing_legend.get_texts()]
        else:
            existing_handles, existing_labels = [], []
        
        scatter1 = addMarker(inset_ax, model, occ, "blue", "spin axis 1", 2)
        
        existing_handles.append(scatter1)
        existing_labels.append("spin axis 1")
        
        if not isNone(model2):
            scatter2 = addMarker(inset_ax, model2, occ, "m", "spin axis 2", 2)
            
            existing_handles.append(scatter2)
            existing_labels.append("spin axis 2")
    
        inset_ax.set_aspect('equal')
        #ax.legend(existing_handles, existing_labels, loc = "upper right")
    
    def plot_fit_cords(self, occ, projected_vertices, edges, title = "", ax = None, block = False, centroide = None):
  
        if ax == None:
            fig, ax = plt.subplots(figsize=(8, 8))
    
        # Tracer les cordes observées
        self.plot_cords_with_errors(occ, "red", ax, centroide = centroide)
    
        # Tracer la silhouette ajustée
        self.plot_silhouette(projected_vertices, edges, "blue", ax, centroide = centroide)
    
        ax.set_xlabel('η (km)', fontsize = 24)
        ax.set_ylabel('ξ (km)', fontsize = 24)
        if title == "":
            ax.set_title('Silhouette ajustée superposée aux cordes observées')
        else:
            ax.set_title(title)
       
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.grid(True)
        ax.set_aspect('equal', adjustable='box')
        plt.show(block=block)  
        
        return ax
            
        
    
    def plot_cords_with_errors(self, occ, c='red', ax=None, block=False, centroide = None):
        
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))
    
        if isinstance(centroide, type(None)):
            centroidex = 0
            centroidey = 0  
            centroide = np.zeros((2))
        else:
            centroidex = centroide[1]
            centroidey = centroide[0]
            centroide = np.asarray(centroide)
            
        # Boucle pour tracer les cordes avec leurs barres d'erreurs
        for i, cordId in enumerate(occ.cords):
            cord = occ.cords[cordId]  # Positions principales de la corde
            status = occ.cordsStatue[cordId]
            cord_err = occ.cordsErr[cordId]  # Erreurs pré-calculées
    
            if status == 0:  # Observation positive
                y_values = cord[:, 0]
                x_values = cord[:, 1]
                
                y_values = y_values - centroidey
                x_values = x_values - centroidex
                
                # Tracer la corde principale
                #ax.plot(x_values, y_values, 'r-', label='Chords' if i == 0 else "")
                ax.plot(x_values, y_values, 'r-')
                ax.scatter(x_values, y_values, marker="x")
                ax.text(np.mean(x_values), np.mean(y_values), f'{cordId}', color=c)
    
                # Tracer les barres d'erreurs en utilisant les positions stockées dans cordsErr
                # Disparition (D)
                d_err_p = cord_err['d_err_p'] - centroide
                d_err_n = cord_err['d_err_n'] - centroide
                # Tracer une ligne entre d_err_p et d_err_n pour représenter l'incertitude
                ax.plot([d_err_p[1], d_err_n[1]], [d_err_p[0], d_err_n[0]], color='blue', linestyle='-', linewidth=1)
    
                # Réapparition (R)
                r_err_p = cord_err['r_err_p'] - centroide
                r_err_n = cord_err['r_err_n'] - centroide
                # Tracer une ligne entre r_err_p et r_err_n pour représenter l'incertitude
                ax.plot([r_err_p[1], r_err_n[1]], [r_err_p[0], r_err_n[0]], color='blue', linestyle='-', linewidth=1)
    
    
            elif status == 1:  # Observation négative
                y_values = cord[:, 0]
                x_values = cord[:, 1]
                
                y_values = y_values - centroidey
                x_values = x_values - centroidex
                
                mean_dir = occ.mean_direction
                slope = mean_dir[0] / mean_dir[1] if mean_dir[1] != 0 else 0
                b = y_values[0] - slope * x_values[0]
                x = [x_values[0] - 15, x_values[0] + 15]
                y = [slope * x[0] + b, slope * x[1] + b]
                #ax.plot(x, y, 'r--', label='negative observations' if i == 0 else "")
                ax.plot(x, y, 'r--')
                ax.scatter(x_values, y_values, marker="x")
                ax.text(x_values[0], y_values[0], f'{cordId}', color=c)
                
            elif status == 2: #Observation ou la disparition ou la réaparition n'a pas été observé
              
                if isinstance(cord[0, 0], type(None)):
                    y = cord[1, 0] - centroidey
                    x = cord[1, 1] - centroidex
                    err_n = cord_err["r_err_n"] - centroide
                    err_p = cord_err["r_err_p"] - centroide
                     
                else:
                    y = cord[0, 0] - centroidey
                    x = cord[0, 1] - centroidex
                    err_n = cord_err["d_err_n"] - centroide
                    err_p = cord_err["d_err_p"] - centroide
                
                print(cord_err)
                
                ax.scatter(x, y, marker="x")
                ax.plot([err_p[1], err_n[1]], [err_p[0], err_n[0]], color='blue', linestyle='-', linewidth=1)
                ax.text(x, y, f'{cordId}', color=c)
    
        # Configuration du graphique
        ax.set_xlabel('Coordonnée η (km)')
        ax.set_ylabel('Coordonnée ξ (km)')
        ax.set_title('Projection des cordes avec barres d\'erreurs dans le plan fondamental')
        ax.grid(True)
        ax.set_aspect('equal', adjustable='box')
        plt.show(block=block)

        # Set axis labels and aspect ratio
        ax.set_xlabel('Coordonnée η (km)')
        ax.set_ylabel('Coordonnée ξ (km)')
        ax.set_title('Projection des cordes dans le plan fondamental (avec barres d\'erreurs)')
        ax.grid(True)
        ax.set_aspect('equal', adjustable='box')
        plt.show(block=block)
    
    
    def plot_3d_chart(self, occ, shouldPlotSolutionOne = True):
        
        
        
       def set_equal_aspect(ax):
            
            extents = np.array([getattr(ax, f'get_{dim}lim')() for dim in 'xyz'])
            sz = extents[:,1] - extents[:,0]
            centers = np.mean(extents, axis=1)
            maxsize = max(abs(sz))
         
            for ctr, dim in zip(centers, 'xyz'):
                getattr(ax, f'set_{dim}lim')(ctr - maxsize/2, ctr + maxsize/2)
       
       t = occ.t_ref()
      
       if shouldPlotSolutionOne:
           self.__setAngle(t, self.model1)
           model = self.model1
       else:
           self.__setAngle(t, self.model2)
           model = self.model2
           
           
       tris = model.triangles
       solution = model.solution

       model.rotate_along_spin(solution["angle"], self.shouldPorject)
       #verts, _, _ = self.compute_projection(solution["lambda"], solution["beta"], solution["angle"], model, occ)
       verts, _, _ = self.compute_projection2(solution["lambda"], solution["beta"], t, solution["period"], model, occ)
        
       fig, ax = occ.plot3dChart(False, True)
       # Exclure les 10 premiers triangles du modèle gris
       num_triangles_to_exclude = min(10, len(tris))
       tris_grey = tris[num_triangles_to_exclude:]  # Triangles à dessiner en gris (à partir du 11e)
       tris_blue = tris[:num_triangles_to_exclude]  # Les 10 premiers triangles à dessiner en bleu
    
       ax.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                       triangles=tris_grey, color='grey', edgecolor='none', alpha=0.7)
    
       ax.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                       triangles=tris_blue, color='blue', edgecolor='none', alpha=0.9)
       
       alpha, delta = self.__compute_geocentric_projection(solution["lambda"], solution['beta'], occ.event.date)
       print("alpha", np.rad2deg(alpha), np.rad2deg(delta))
       xyz = self.__spherical2cartesian(alpha, delta)
            
       ax.quiver(0, 0, 0, xyz[0], xyz[1], xyz[2], color='red', label='Axe de rotation', hatch="-")  # Axe de rotation

       ax.axis('off')
       ax.legend()
       ax.set_facecolor('black')
       

       
       plt.gca().set_aspect('equal', adjustable='box')
       
       fig.patch.set_facecolor('black')
       plt.show(block=False)

    

def pathFormat(p):
    
    if p[0] == ".":
        p = file_path + p[1::]
        
    return p

def argument():
    parser = argparse.ArgumentParser(description="""This script offert 4 tools. The first tool is a parser for xml file link to stellar occultations (the class Event is in charge).
                                     The second tool is a model parser usefull to parse file model.tri. this tool is able to produce multiple plot of the asteroid model and can
                                     replace the code 'draw_shape.m'. the classe in shape is ModelParser. The third tool is manage by the class Occultation. this tool is highly related to 
                                     the first tool as it's take as intupe an object Event and allow to plot cords for exemple. The last tool is managed by the class OccultationFit and is 
                                     able to fit an occultation according to a model. if you want to use this code as script, please read the following description and options. 
                                     If you want to add this project to your code and import it in your script each classes and tools of this code are described here: https://github.com/antoinech13.

                                     These script ofer two type of options. the frist type describe the tool that you want to use. second type of options are options link to the tool that you are using.
                                     some of these options are the same between tools which can add ambiguity if you decide to use multiple tools in one commande. for these reason, the usage of only one tool is allow.
                                     More that one tool will raise an error and the code will stop. options to use tool will start by 'TOOL' in the description and all other options which can be use will also be writen
                                     so please pay attention to options:""")

    parser.add_argument("-m", "--modelTool",  action = "store_true", help= "TOOL: this option is dedicated to plot parse and plot model.tri files. should be use with -mP and ")
    parser.add_argument("-oF", "--occultationFit", action = 'store_true', help = "TOOL: this option is dedicated to fit and plot occultation with model")
    parser.add_argument("-o", "--occultation", action = 'store_true', help = "TOOL: this option is dedicated to build cord from events.")
    
    parser.add_argument("-f", "--fit", action = "store_true", help = "If use if -oF: launch the fitting algorithm")
    
    parser.add_argument('-p', "--plotModel", action = "store_true", help = """if use with -m: plot the 3d model
                                                                              if use with -o: display the observation
                                                                              if use with -of: display geometry of the occultation""")
                                                                              
    parser.add_argument('-drawShape', "--drawShape", action = "store_true", help = "plot 3 figures of the 3d model. do the same than draw_shape.m")
    parser.add_argument("-mP", "--modelPath", type=str, default = file_path, help = """IF USE WITH -m : specify the path of model.tri files. if path of model.tri file use, the code will read this specific file. if a directory
                                                                                       path is use, the code will search for all model.tri file in this directory, DEFAULT = the directory of the script
                                                                                       IF USE WITH -oF : specify the path of model.tri files for the first solution. If only one solution is availabe, this option should be use. 
                                                                                       IF NONE OF -mP OR -m2P ARE USED THE CODE WILL SEARCH FOR *1.tri AND *2.tri""")
    parser.add_argument("-m2P", "--model2Path", type=str, default = file_path, help = "USE WITH -oF: specify the path of model.tri files for the second solution." )
    
    parser.add_argument("-input", "--input", type=str, default = file_path + "ocFit_input", help = "path to input file of occultation fit")
    parser.add_argument("-ocP", "--occultationPath", type=str, default = file_path, help = "directory where to find occultations. Can be with with -o and -oF")
    
    parser.add_argument("-sD", "--shouldRequestStarFromDatabase", action = "store_true", help = "If this flag is present, the script will ignore stars cordinnates written in XML occultations files and request it from database. Usable with -o and -oF")
    parser.add_argument("-sI", "--shouldBeIndependanteFiting", action = "store_true", help = "If this flag is present, the script will fit each occultation independently")
    parser.add_argument("-sT", "--solverType", type=str, default = "SLSQP", help = "Define the solver type from scipy. default = SLSQP. ")
    parser.add_argument("-sS", "--solverStartegy", type=str, default = "hybrid", help = """Define the solver strategy. 6 strategies exist (both deep strategies are not used for uncertainties due to time of execution. will be replace by the corresponding hybrid technic) :  - scipy (use only scipy solver) \n
                                                                                                                                                                                                                                                                                 - home (use homemade solver) \n
                                                                                                                                                                                                                                                                                 - hybrid (use scipy first then homemade) \n
                                                                                                                                                                                                                                                                                 - deep_hybrid (use scipy first then run multiple time homemade solver with small different initial point from the first run result) \n 
                                                                                                                                                                                                                                                                                 - hybrid_r (use first homemade solver than scipy. reverse of hybrid) \n
                                                                                                                                                                                                                                                                                 - deep_hybrid_r (reverse of deep_hybrid)\n""") 
    parser.add_argument("-nSd", "--number_of_sub_div", type=int, default = 2, help = """"Define the number of subdivision that the home optimizer do. A subdivision is a full optimisation 
                                                                                         process until reaching an optimum. for each subdivision, the steps used are divided by the subdivision index. 
                                                                                         For example, if this parameter is set to 3 than the optimizer will optimize a first time with a variable step 
                                                                                         until not moving or until triggering the -mSs than it will optimise it one more time with step/10 and a last time with step/100.
                                                                                         higher number mean more accurate results but also more computational time. DEFAULT= 2""")
    parser.add_argument("-mSs", "--max_sub_stop", type = int, default = 75, help = """Fixe a limit of duration for a subdivision. A subdivision will iteratively optimise the code until 
                                                                                      finding a fix optimum. If they are not fix optimum after -mSs iteration, than the subdivision stop anyway. DEFAULT = 75""")
    parser.add_argument("-noErr", "--no_err", action = 'store_true', help = "If use with -oF, ignore cords Errors for fit")
    parser.add_argument("-pC", "--plotchord", action = 'store_true', help = "If use with -o, display cords")
    parser.add_argument("-pObs", "--plotObservation", action = 'store_true', help = "If use with -o, display the observation")
    parser.add_argument("-pS", "--plotSilhouette", type = float, default = None, help = """Take an float as input (angle or rotation of the asteroid along it axis) 
                                                                                           If use with -oF, display the silhouette obtain at the given angle.""")
    
    parser.add_argument("-pI", "--projectAlongInertia", action = "store_true", help = """If use with -oF: fitting of occultation usually consider that the z axis of model is the axis of rotation which is usually true. 
                                                                                         But in some case the main axis of rotation (computed from tensior of inertia) is not exactly aline with the z axis.
                                                                                         The difference is really small but in case where this difference matter, use this flag to projet vertices of model
                                                                                         along eigenvectors of the tensor of inertia
                                                                                         If use with -m: eigen vectors are usually not unique. Therefore, nativly, the parameter -drawShape will organise them to 
                                                                                         obtain the same result than the drawshape.m script which is not all time the best way (mirror effect, axis inversed, etc...). This parameter will organise eigenvectors
                                                                                         to ensure that main axis are all positif and in the diagonal of the eigenvectors matrix""")
                                                        
    parser.add_argument("-backEnd", "--backEnd", type = str, default = "qt", help = "change graphic back end for matplotlib. if any error are trigger you can change this parameter options are agg, tkagg, qt. DEFAULT: qt")
    
    args = parser.parse_args()                                                                               
    
    return args
    
    


def modelTool(args):
    
    shouldPlotMod = args.plotModel 
    shouldDrawShape = args.drawShape   
    modelPath = args.modelPath
    

    if "*.tri" in modelPath:
        modelPath = [modelPath]
    else:
        modelPath = glob.glob(modelPath + "*.tri")
    
    for file in modelPath:
        print("model of file: ", file)
        mp = ModelParser(file)
        
        if shouldPlotMod:
            mp.plot_vertices(block = True)
        
        if shouldDrawShape:
            mp.plot3_vertices(asMatlab = not args.projectAlongInertia, block = False)
            mp.plot3_vertices(True, not args.projectAlongInertia, False)
            input("press a touch to check next .tri file")
            plt.close()
    
    
    input("press touch to exit")
    
    sys.exit()

def isNone(val):
    return isinstance(val, type(None))

def fitOc(args):
    
    def getModelPath(p, solution):
        
        p = pathFormat(p)
        
        if os.path.isfile(p):
            return p
        
        elif os.path.isdir(p):
            p = glob.glob(p + "*"+ str(solution) + ".tri")
            
            if len(p) == 0:
                print("Warning: no model *.tri was found for solution", solution)
                return None
            
            p = p[0]
            
            print("WARNING ! no model.tri was provide," , p, "was found and will be use as solution", solution)
        
        else:
            print("ERROR: no model *.tri was found for solution", solution)
            sys.exit()
            
        return p
    
    pa_input = pathFormat(args.input)
    ofi = OcFitInputParser(pa_input)
    
    m1 = None
    m2 = None

    m1_path = getModelPath(args.modelPath, 1) 
    m2_path = getModelPath(args.model2Path, 2)
    #m1_path = "./aegle/Aegle/flat/s1/0.005/model.tri"
    #m2_path = "./aegle/Aegle/flat/s2/0.0003/model.tri"
    
    if isNone(m1_path) and isNone(m2_path):
        print("ERROR: no .tri file available.")
        sys.exit()

    if ofi.nSolution == 2 and ( isNone(m1_path) or isNone(m2_path) ):
        print("ERROR: two solutions in input file but only one model provided")
        sys.exit()
    
    if ofi.nSolution == 1 and not isNone(m1_path) and not isNone(m2_path):
        print("Warning ! one solution provided in inputfile but two models found. the first model will be use.")
        
    if isNone(m1_path) and not isNone(m2_path):
        m1_path = m2_path
        m2_path = None
        
    s1 = {"lambda" : ofi.getLamb(0),
          "beta": ofi.getBeta(0),
          "period": ofi.getPeriod(0),
          "zero": ofi.getInitialTime(0),
          "phase": ofi.getInitialPhase(0),
          "err_l_m" : ofi.getLambErr(0)[0],
          "err_l_p" : ofi.getLambErr(0)[1],
          "err_b_m": ofi.getBetaErr(0)[0],
          "err_b_p": ofi.getBetaErr(0)[1],
          "err_p_m" : ofi.getPeriodErr(0)[0],
          "err_p_p" : ofi.getPeriodErr(0)[1],
          "shouldLambdaBeOptimised" : ofi.l_opt[0],
          "shouldBetaBeOptimised" : ofi.b_opt[0],
          "shouldPeriodBeOptimised" : ofi.p_opt[0]
          }
    
    
    m1 = ModelParser(m1_path, s1)
     
    if ofi.nSolution == 2:
        
        s2 = {"lambda" : ofi.getLamb(1),
              "beta": ofi.getBeta(1),
              "period": ofi.getPeriod(1),
              "zero": ofi.getInitialTime(1),
              "phase": ofi.getInitialPhase(1),
              "err_l_m" : ofi.getLambErr(1)[0],
              "err_l_p" : ofi.getLambErr(1)[1],
              "err_b_m": ofi.getBetaErr(1)[0],
              "err_b_p": ofi.getBetaErr(1)[1],
              "err_p_m" : ofi.getPeriodErr(1)[0],
              "err_p_p" : ofi.getPeriodErr(1)[1],
              "shouldLambdaBeOptimised" : ofi.l_opt[1],
              "shouldBetaBeOptimised" : ofi.b_opt[1],
              "shouldPeriodBeOptimised" : ofi.p_opt[1]
              }
    
        m2 = ModelParser(m2_path, s2)
            
    evs = []
    occs = []
    
    paths = glob.glob(args.occultationPath + "*.xml")
    
    if len(paths) == 0:
        print("No xml files found")
        sys.exit()
    
    for path in paths:
        evs.append(Event(path, args.shouldRequestStarFromDatabase))
    
    #=============correct shift==============================
 
    for key in ofi.cordToShift:
        for ev in evs:
            if ev.date.to_value("iso").split(" ")[0] == key:
                for element in ofi.cordToShift[key]:
                    if not isinstance( ev.observer[element[0]]["d"], type(None)):
                        ev.observer[element[0]]["d"] = ev.observer[element[0]]["d"] + element[1] / (24 * 3600)
                    if not isinstance( ev.observer[element[0]]["r"], type(None)):
                        ev.observer[element[0]]["r"] = ev.observer[element[0]]["r"] + element[1] / (24 * 3600)
    #========================================================
    for ev in evs:
        occs.append(Occultation(ev))


    #===========delete bad chords===========================
  
    for key in ofi.cordToDelete:
        for occ in occs:
                if occ.event.date.to_value("iso").split(" ")[0] == key:
                    for element in ofi.cordToDelete[key]:
                        del occ.cords[element]
                        del occ.cordsStatue[element]
                
    #========================================================
   
    #===========cord to weight===============================
   
    for key in ofi.weight:
        for occ in occs:
            if occ.event.date.to_value("iso").split(" ")[0] == key:
                for element in ofi.weight[key]:
                    occ.cordWeight[element[0]] = element[1]
                
    
    
    
    ocf = OccultationFit(occs, m1, m2, args.projectAlongInertia)    
    
    if args.fit:
        ocf.fit(args.shouldBeIndependanteFiting, args.solverType, args.solverStartegy, args.number_of_sub_div, args.max_sub_stop, ofi.multi, ofi.centroid, ofi.size, args.no_err)
        
    if args.plotModel:
        for occ in occs:
            ocf.plot_3d_chart(occ)
            if not isNone(m2):
                ocf.plot_3d_chart(occ, False)
   
    if not isNone(args.plotSilhouette):
        
        ang = args.plotSilhouette
        lamb = m1.solution["lambda"]
        beta = m1.solution["beta"]
        period = m1.solution["period"]
        print(lamb, beta, period)
        for occ in occs:
            m1.rotate_along_spin(ang)
            v, pv, e = ocf.compute_projection2(lamb, beta, occ.t_ref(), period, m1, occ)
            ocf.plot_silhouette(pv,e )
        
    input("press touch to exit")
    sys.exit()
    
def occultation(args):
    
    paths = glob.glob(args.occultationPath + "*.xml")
    
    for path in paths:
        occ = Occultation(Event(path, args.shouldRequestStarFromDatabase))
        
        if args.plotchord:
            occ.plot()
            
        if args.plotModel:
            occ.plot3dChart()
    
    input("press touch to exit")
    sys.exit()
    
def parsArguments(args):
    
    
    if args.backEnd.lower() == "qt":       
        print("qt chose")
        matplotlib.use('Qt5Agg')
    elif args.backEnd.lower() == "agg":
        print('agg chose')
        matplotlib.use('Agg')
    elif args.backEnd.lower() == "tkagg":
        print("tkadd chose")
        matplotlib.use('TkAgg')
    
    if args.modelTool:
        modelTool(args)
        
    elif args.occultationFit:
        fitOc(args)
    
    elif args.occultation:
        occultation(args)
        
if __name__ == "__main__":
    
    
    args = argument()
    parsArguments(args)
    