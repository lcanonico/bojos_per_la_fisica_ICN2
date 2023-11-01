import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import kwant
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from numpy.random import random
welcome="Bibliotecas cargadas. ¡Bienvenidos!"

import requests
import IPython.display as Disp
grafeno_real = 'https://i.pinimg.com/564x/27/73/60/277360f8a8265b0a8008d1d03023a7be.jpg'
grafeno_disp = 'https://icn2.cat/images/20190802_spin-communication.jpg'
def ejemplo_grafeno_real():
  return Disp.Image(requests.get(grafeno_real).content);
def ejemplo_dispositivo_grafeno():
  return Disp.Image(requests.get(grafeno_disp).content);


#creacion del grafeno
a0 = 0.246 #[nm]
acc = a0/np.sqrt(3) #[nm]
latvecs = [(a0,0),(0,3*acc)]
latpos=[(0,0),(0,acc),(0.5*a0,3*0.5*acc),(0.5*a0,5*0.5*acc)]
graphene = kwant.lattice.general(latvecs,latpos,norbs=1)
suba1,subb1,suba2,subb2 = graphene.sublattices
t=-2.71
energy = -2;


def create_system(L, W, m=0,U=1.0,c=0.0,sym=None,r0=(0,0)):
  syst = kwant.Builder()
  if sym is not None:
    syst = kwant.Builder(sym)
  def shape(pos):
    x,y = pos
    return (-L/2<= x <=L/2) and (-W/2<=y<=W/2)
  def onsite(site):
    if random()*100<=c:
      if(site.family == suba1 or site.family==suba2):
        return m + U*(random()-0.5)
      else:
        return -m + U*(random()-0.5)
    else:
      if(site.family == suba1 or site.family==suba2):
        return m 
      else:
        return -m
  syst[graphene.shape(shape,r0)] = onsite
  ##Within the same cell
  syst[kwant.builder.HoppingKind((0,0),subb1,suba1)]=t
  syst[kwant.builder.HoppingKind((0,0),subb2,suba2)]=t
  syst[kwant.builder.HoppingKind((0,0),subb1,suba2)]=t
  #Fuera de la celula unitaria
  syst[kwant.builder.HoppingKind((1,0),subb1,suba2)]=t
  syst[kwant.builder.HoppingKind((0,-1),subb2,suba1)]=t
  syst[kwant.builder.HoppingKind((-1,-1),subb2,suba1)]=t
  syst.eradicate_dangling()
  return syst
  

    




"""energy = -2;
lat_c  = 0.246;
lat_vec= lat_c*np.array(((1, 0), (0.5,0.5*np.sqrt(3))));
Area = np.sqrt(3)*0.5*lat_c*lat_c;
orbs   = lat_c*np.array([(0, 0), (0, 1 / np.sqrt(3))]);
graphene = kwant.lattice.general(lat_vec, orbs);
a, b = graphene.sublattices

def create_system( L, W, sym=None, U=1.0, c=0.0, r0=(0,0), phi=0 ):

  syst = kwant.Builder()
  if sym is not None:
    syst = kwant.Builder(sym);
    phi = 0
    
  def shape(pos):
    x, y = pos;
    a0, a1 = np.linalg.inv(lat_vec).T.dot(pos);
    return ( r0[0]/lat_c <= a0 <= L/lat_c ) and ( r0[1]/lat_c <= a1 <= W/lat_c )

  def anderson(site):
      if random()*100 <= c:
        return U*(random()-0.5);
      return 0.0;
  #incorporate anderson disorder as onsites
  syst[graphene.shape(shape, r0)] = anderson;

  def hopping(site_i, site_j):
    xi, yi = site_i.pos;
    xj, yj = site_j.pos;
    return -2.8*np.exp(-0.5j * phi * (xi - xj) * (yi + yj))
  #incorporate hoppings
  syst[graphene.neighbors()] = hopping;

  return syst;"""

def crear_cable(L, W,m=0, U=1.0, c=0.0):
  return create_system(L, W, m,U,c,sym=None)
  
def crear_cable_infinito(L, W,m=0, U=1.0, c=0.0):
  dir = latvecs[0]
  s = kwant.TranslationalSymmetry(dir)
  return create_system(L, W, m,U,c,sym=s)

def agregar_contactos(syst,L, W):
    tdir=-graphene.vec((1,0));
    sym = kwant.TranslationalSymmetry(tdir);
    lead = create_system(L=L, W=W, sym=sym,r0=2.0*tdir);
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())
    return syst;

def family_color(site):
    if site.family == suba1 or site.family==suba2:
        return 'red'
    else: 
        return 'blue'
    

def graficar_sistema( syst, ax=None ):
  if(ax):
    kwant.plot(syst, show=False,dpi=100,site_color=family_color,site_lw=0.1,ax=ax);
  else:
    kwant.plot(syst, show=False,dpi=100,site_color=family_color,site_lw=0.1);
  imp_pos = [s.pos for s,v in syst.site_value_pairs() if np.abs(v(s))>1e-5 ]
  if len(imp_pos)>0: 
    ax = plt.gca();
    ax.scatter( *np.transpose(imp_pos), c="k",zorder=3,s=1);
  plt.show();

def grafica_densidad(syst,E0):
  fsyst = syst.finalized()
  wf = kwant.wave_function(fsyst,energy=E0)
  psi = wf(0)[0]
  rho = kwant.operator.Density(fsyst)
  ldos = rho(psi)
  kwant.plotter.density(fsyst,ldos,cmap=mpl.colormaps["viridis"])
  plt.show()

def calcula_conductancia(fsyst,E0,nreal=1 ):
  C0 = 7.7480e-5; #conductance quantum
  return C0*np.mean([ kwant.smatrix(fsyst, E0).transmission(1, 0) for i in range(nreal)]); 

def calcula_matriz_conductancia(fsyst,E0,nreal=1 ):
  C0 = 7.7480e-5;#conductance quantum
  C = C0*kwant.smatrix(fsyst,E0 ).conductance_matrix();
  return C[:3,:3]; 

def calcula_matriz_resistencia(fsyst,E0,nreal=1 ):
  return np.linalg.inv(calcula_matriz_conductancia(fsyst,E0,nreal)); 


def calcula_resistencia(fsyst,E0,nreal=1 ):
  return 1/calcula_conductancia(fsyst,E0,nreal );

def calcula_corriente(syst,sample=50,nreal=1):
  C0 = 7.7480e-5#conductance quantum
  fsyst = syst.finalized()
  energies = np.linspace(-1,1,sample)
  conductance = np.zeros(sample,dtype=float)
  for r in range(0,nreal):
    for idx,e in enumerate(energies):
        smatrix = kwant.smatrix(fsyst,e)
        conductance[idx]+=C0*smatrix.transmission(1,0)/nreal
  return energies, conductance


def calcula_bandas(syst,sample=50):
  fsyst = syst.finalized()
  ebands = kwant.physics.Bands(fsyst)
  momenta = np.linspace(0, np.pi*2.,sample)
  energies = []
  for k in momenta:
    energies.append(ebands(k))
  return momenta,energies



