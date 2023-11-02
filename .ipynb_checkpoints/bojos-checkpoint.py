import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import kwant
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection

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


def intersectar_circulos(r1=1, r2=1, d=0):
    assert r1 >= 0 and r2 >= 0 and d >= 0, r'introduzca una distancia positiva'
    r1, r2 = np.max([r1, r2]), np.min([r1, r2])
    
    if d >= r1 + r2:
        return 0
    elif d <= r1 - r2:
        return np.pi * r2 ** 2
    else:
        d1 = (r1 ** 2 - r2 ** 2 + d ** 2) / (2 * d)
        d2 = d - d1
        return r1 ** 2 * np.arccos(d1 / r1) - d1 * np.sqrt(r1 ** 2 - d1 ** 2) + \
               r2 ** 2 * np.arccos(d2 / r2) - d2 * np.sqrt(r2 ** 2 - d2 ** 2)

def hamiltoniano_nanocable1D(k, e1=0, e2=0, t11=0, t22=0, t12=0):
    return np.array([[e1+2*t11*np.cos(k), t12*np.exp(1j*k)], [t12*np.exp(-1j*k), e2+2*t22*np.cos(k)]])

def graficar_acoplamientos_nanocable1D(radio1, radio2, distancia, e1=-10, e2=0):
    t11, t12, t22 = intersectar_circulos(r1=radio1, r2=radio1, d=distancia), intersectar_circulos(r1=radio1, r2=radio2, d=distancia), intersectar_circulos(r1=radio2, r2=radio2, d=distancia)

    fig = plt.figure(figsize=(10, 2.5))

    ax_orb = fig.add_axes([0, 0, 0.65, 1])
    x_lim = 10
    X = list([-distancia*n for n in range(1, int(x_lim/distancia)+5)[-1::-1]]) + [0] + list([distancia*n for n in range(1, int(x_lim/distancia)+5)])
    Phi = np.linspace(0, 1, 101) * np.pi
    for x in X:
        ax_orb.plot([x], [0], markersize=10, color='k', linestyle='None', marker='o', zorder=0)
        ax_orb.fill_between(x+radio1*np.cos(Phi), y1=radio1*np.sin(Phi), y2=-radio1*np.sin(Phi), color=mpl.colormaps['coolwarm'](0.), alpha=0.5, zorder=-10)
        ax_orb.fill_between(x+radio2*np.cos(Phi), y1=radio2*np.sin(Phi), y2=-radio2*np.sin(Phi), color=mpl.colormaps['coolwarm'](1.), alpha=0.2, zorder=-20)
    ax_orb.set_xlim(-x_lim, x_lim)
    ax_orb.set_xlabel(r'posicion')
    ax_orb.set_title(r'acoplamiento de orbitales atomicos')
    text = f'acoplamiento azul - azul : {t11}\n' + \
           f'acoplamiento rojo - rojo : {t22}\n' + \
           f'acoplamiento azul - rojo : {t12}'
    print(text)
    
    ax_bandas = fig.add_axes([0.75, 0, 0.25, 1])
    K = np.linspace(-1, 1, 501) * np.pi
    Evals, Evecs = np.linalg.eigh(np.array([hamiltoniano_nanocable1D(k, e1=e1, e2=e2, t11=t11, t12=t12, t22=t22) for k in K]), UPLO='U')
    for n in range(2):
        points = np.array([K, Evals[:, n]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, colors=mpl.colormaps['coolwarm'](np.square(np.abs(Evecs[:, 1, n]))), linewidth=3)
        ax_bandas.add_collection(lc)
    ax_bandas.set_xlim(-np.pi, np.pi)
    ax_bandas.set_xticks([-np.pi, 0, np.pi])
    ax_bandas.set_xticklabels(r'$-\pi$ $0$ $+\pi$'.split())
    ax_bandas.set_xlabel(r'$k=\frac{2\pi}{Longitud}$')
    ax_bandas.set_ylim(Evals.min()-1, Evals.max()+1)
    ax_bandas.set_ylabel(r'energia')
    ax_bandas.set_title(r'bandas de energia')
    return fig, [ax_orb, ax_bandas]

def graficar_ondas_nanocable1D(dimerizacion=0, longitud_de_onda=4, N_celdas=5, banda=1):
    assert banda == 0 or banda == 1, 'inserte banda=0 o banda=1'
    
    xA = np.arange(-N_celdas, N_celdas+1)
    xB = np.arange(-N_celdas-1, N_celdas+1) + 0.5 + 0.5 * dimerizacion / 100
    k = 4 * np.pi / longitud_de_onda
    u, v = np.exp(-dimerizacion / 100), np.exp(dimerizacion / 100)
    
    fig = plt.figure(figsize=(10, 2.5))
    
    ax_bandas = fig.add_axes([0, 0, 0.25, 1])
    K = np.linspace(-1, 1, 101) * np.pi
    ax_bandas.plot(K, -np.sqrt(u ** 2 + v ** 2 + 2 * u * v * np.cos(K)), color='k', linewidth=3)
    ax_bandas.plot(K, np.sqrt(u ** 2 + v ** 2 + 2 * u * v * np.cos(K)), color='k', linewidth=3)
    ax_bandas.plot([k], (2*banda-1) * np.sqrt(u ** 2 + v ** 2 + 2 * u * v * np.cos(k)), linestyle='None', marker='o', markersize=10, color='k', label='onda graficada')
    ax_bandas.legend()
    ax_bandas.set_xlim(-np.pi, np.pi)
    ax_bandas.set_xticks([-np.pi, 0, np.pi])
    ax_bandas.set_xticklabels([r'$-\pi$', r'$0$', r'$+\pi$'])
    ax_bandas.set_xlabel(r'$k = \frac{4\pi}{Longitud}$')
    ax_bandas.set_ylabel(r'energia')
    ax_bandas.set_title(r'bandas de energía')
    
    ax_onda = fig.add_axes([0.35, 0, 0.65, 1])
    X = np.linspace(-N_celdas - 0.5, N_celdas + 0.5, 501)
    ax_onda.plot(xA, np.zeros(len(xA)), marker='o', linestyle='None', markersize=5, color='tab:blue')
    ax_onda.plot(xB, np.zeros(len(xB)), marker='o', linestyle='None', markersize=5, color='tab:red')
    [ax_onda.plot([x, x], [0, np.cos(k * (x - 0.5*banda))], marker='None', linestyle='solid', linewidth=3, color='tab:blue') for x in xA];
    [ax_onda.plot([x, x], [0, (1-2*banda) * np.cos(k * (x - 0.5*banda))], marker='None', linestyle='solid', linewidth=3, color='tab:red') for x in xB];
    ax_onda.fill_between(X, y1=np.cos(k * (X - 0.5*banda)), y2=-np.cos(k * (X - 0.5*banda)), linewidth=3, edgecolor='k', facecolor='tab:gray', alpha=0.3, zorder=-100)
    ax_onda.set_xlim(X.min(), X.max())
    ax_onda.set_xticks(np.arange(-N_celdas, N_celdas+1, 0.5))
    ax_onda.set_xticklabels(map(str, (2*np.arange(-N_celdas, N_celdas+1, 0.5)).astype(int)))
    ax_onda.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax_onda.set_yticklabels(r'$-1.0$ $-0.5$ $0.0$ $+0.5$ $+1.0$'.split())
    ax_onda.set_ylabel(r'amplitud')
    ax_onda.set_xlabel(r'posicion')
    ax_onda.set_title(r'visualizacion de una onda en espacio real')
    
    return fig, [ax_bandas, ax_onda]
