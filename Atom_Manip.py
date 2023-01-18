import ReadLAMMPS as rl
import atom
import numpy as np
import os, sys
import pdb
import math
import Plotting as _plot



NA = 6.022e23
Masses = {'1':28.0855,'2': 15.9994,'3': 1.008}

def region(atom):
    d=5.43
    reg=np.empty([3,2])
    for i in [0,1,2]:
        low, high = atom.coords[i] - d, atom.coords[i] + d
        if(low < atom.box[i][0]):
            low=atom.box[i][0]
        elif(high > atom.box[i][1]):
            high=atom.box[i][1]
        reg[i] = [low,high]
    return reg


def get_env(at):
    i=0
    nlist=[]
    simask=[]
    for n in at.neighbors:
        nlist.append(atom.distance(at,n))
        simask.append(n.type==2)
    nlist=np.asarray(nlist) #array of at neighbor distances
    simask=np.asarray(simask) #Bool array of at silicon neighbors
    bmask=np.logical_and(nlist>1.3,nlist<1.6) #bool array for is at bonded to a Si
    nnmask=np.logical_and(nlist>2.0, nlist<2.8) #bool array for how many neighbors are 2.4->2.8 A away
    sinn=np.logical_and(simask,nnmask) #bool array for are Si neighbors between 2.4->2.8 
    bsi=np.logical_and(simask,bmask) #bool array for is the neighbor 1.4->1.6 a Si atom
    SiBonded=(np.count_nonzero(bsi)==1) #count to make sure only one Si neighbor is 1.4->1.6, if more than 1
    # then H atom is at a bond center.
    count_nn = np.count_nonzero(sinn)
    if(SiBonded and count_nn >= 1):
        nn = at.neighbors[sinn]
        return move_H(at, nn, count_nn)

    else:
        return np.asarray([at])

def move_H(at,nn, cnt):
    if(cnt >= 2):
        for n in nn:
            for oth in nn:
                if(n.id != oth.id):
                    dist=atom.distance(n,oth)
                    if(dist > 2.15 and dist < 2.8):
                        return move_bc(at,n,oth)
    #if no nn are also neighbors then move H to an interstitial site away from the two si atoms.
    return move_inter(at,nn)

def move_bc(at, n1, n2):
    print("BC")
    r=n1.coords - n2.coords
    r_hat=r/np.sqrt(np.sum(r*r))
    mid=(n1.coords + n2.coords)/2.0
    at.coords=mid
    n1.coords= n1.coords + 0.45*r_hat
    n2.coords= n2.coords - 0.45*r_hat
    at.distance_moved()
    n1.distance_moved()
    n2.distance_moved()
    return np.array([at,n1,n2])

def move_inter(at, nn):
    s=np.size(nn)
    if(s >= 2):
        r1= nn[0].coords - at.coords
        r2= nn[1].coords - at.coords
        r_vec= r1 + r2
        r_hat=r_vec/np.sqrt(np.sum(r_vec*r_vec))
        at.coords = at.coords + 3.35*r_hat
        print("2 Neighbors")
        at.distance_moved()
        return np.asarray([at])
    elif(s==1):
        r_vec= nn[0].coords - at.coords
        r_hat=r_vec/np.sqrt(np.sum(r_vec*r_vec))
        at.coords = nn[0].coords + 1.5*r_hat
        print("1 Neighbor")
        at.distance_moved()
        return np.asarray([at])
    else:
        return np.asarray([at])

def get_hydrogen(array_of_atoms):
    mask=np.empty([np.size(array_of_atoms)],dtype=bool)
    i=0
    cont=0
    for at in array_of_atoms:
        if(at.type == 1):
            mask[i]=True
            i+=1
            cont+=1
        elif(at.type == 2):
            mask[i]=False
            i+=1
    return array_of_atoms[mask]

def Migration(at, file):
    print_env(at)
    moved=get_env(at)
    final = io.write_final(moved, at.id, file)
    if not(final):
        print("Hydrogen atom {0} was not moved")
        return 1
    else:
        io.comp(file,final)
        return 0

def print_env(at):
    atype=["","H","Si"]
    print("")
    txt1="{0} {1}: {2:.2f} {3:.2f} {4:.2f}".format(atype[at.type], at.id, *at.coords)
    txt2="{0} {1} {2:.3f}"
    print(txt1)
    print("Type: ID: Distance:")
    for n in at.neighbors:
        ftext=txt2.format(atype[n.type], n.id, atom.distance(at,n))
        print(ftext)
    return 2



def Compute_Density(array_of_atoms):
    mass = []
    cut = 9.8
    vol = atom.Compute_Volume(cut)
    for at in array_of_atoms:
        if(at.coords[2] > cut):
            mass.append(Atom_Masses[at.type])
    num = len(mass)
    M_tot = np.sum(mass)
    Mcm = M_tot/NA
    rho = Mcm/vol
    print("Num = {0:.0f}, M(g)= {1:.2e}, Vol = {2:.2e}, rho = {3:.2f}".format(num,Mcm, vol, rho))
    return num, rho


def H2_List(H_ats):
    lst = []
    for at in H_ats:
        if(np.isin(at,lst)):
            continue
        else:
            for neb in at.neighbors:
                if(neb.type == 1 and atom.distance(at,neb) <= 1.0):
                    lst.append([at,neb])
                else:
                    continue
    return lst

def Four_Coord(at):
    if(len(at.neighbors) >= 4):
        for neb in at.neighbors:
            if(neb.type == 1):
                return False
            else:
                continue
        return True
    else:
        return False

def Create_DB(H2, md):
    out_file = out_path + "/DB-S{0}-H{1:.0f}.dump".format(ARGS[2],H2[0].id)
    center = H2[0].coords + 0.5*atom.Compute_Vec(H2[0], H2[1])
    f_init = [1000,1,0,0,0]
    fake = atom.ATOM(f_init)
    fake.coords = center
    tmp = md.atoms
    at_list = []
    idx = 10000
    min_dist = 4.0
    for at in md.atoms:
        curr_dist = atom.distance(fake,at)
        if(curr_dist < min_dist and Four_Coord(at)):
            min_dist = curr_dist
            idx = at.id
    print("\n###################\n{0:.0f}\n".format(idx))
    md.pop(idx)
    io.Write_Dump(out_file,md)
    md.atoms = tmp
    md.num +=1
    return


def Inspect_Region(sim, region):
    _dict = {'1' : 0 , '2': 0, 'rho': 0}
    vol = 1.0
    mass = 0.0
    for i in range(3):
        vol *= abs(region[i][1] - region[i][0])
    for at in sim.atoms:
        if(at.IS_IN(region)):
            _dict[str(at.type)] +=1
            mass += Masses[str(at.type)]
    _dict['rho'] = mass/(NA*vol*1.0e-24)
    return _dict


def Compute_Ratio(sim):
    #(1) break up the simulation box in to rectangular parallapiped i.e. R = w X w X h for w = 5.43 => 25 regions.
    # what I really want is number of bins x and number of bins y
    #pdb.set_trace()
    num_x = 16
    num_y = 10
    dx = sim.box[0][1]/num_x
    dy = sim.box[1][1]/num_y
    #then loop through atoms checking their x and y values and then assign to a region
    counts_Si = np.full([num_x,num_y], -24)
    counts_Ox = np.zeros([num_x,num_y])
    temp = 0
    for at in sim.atoms:
        x, y, z = at.coords
        x_idx = math.floor(x/dx)
        y_idx = math.floor(y/dy)
        if(x_idx > num_x - 1):
            x_idx = num_x -1
        if(y_idx > num_y -1):
            y_idx = num_y -1
        if(at.type == 1 ):
            counts_Si[x_idx][y_idx]+=1
        elif(at.type == 2):
            counts_Ox[x_idx][y_idx]+=1
    hist = []
    for i in range(num_x):
        hist_x = i*dx + dx/2.0
        for j in range(num_y):
            hist_y = j*dy + dy/2.0
            nsi = counts_Si[i][j]
            nox = counts_Ox[i][j]
            tot = nsi + nox
            rat = nox/nsi
            hist.append([hist_x, hist_y, nsi, nox,rat, tot])
    return np.asarray(hist)



def Compute_Ratio2(sim):
    #(1) break up the simulation box in to rectangular parallapiped i.e. R = w X w X h for w = 5.43 => 25 regions.
    # what I really want is number of bins x and number of bins y
    #pdb.set_trace()
    num_x = 16 #int(sim.box[0][1]/5.43)
    num_y = 10 #int(sim.box[1][1]/5.43)
    dx = sim.box[0][1]/num_x
    dy = sim.box[1][1]/num_y
    #then loop through atoms checking their x and y values and then assign to a region
    counts_Si = np.zeros([num_x,num_y])
    counts_Ox = np.zeros([num_x,num_y])
    temp = 0
    for at in sim.atoms:
        x, y, z = at.coords
        if(z < 16.0):
            continue
        idx = math.floor(x/dx)
        idy = math.floor(y/dy)
        test1 = round(x/dx)
        test2 = round(y/dy)
        #pdb.set_trace()
        if(idx > num_x - 1):
            idx = num_x -1
        if(idy > num_y -1):
            idy = num_y -1
        if(at.type == 1 ):
            counts_Si[idx][idy]+=1
        elif(at.type == 2):
            counts_Ox[idx][idy]+=1

    ratio = np.zeros([num_x,num_y])
    
    for i in range(num_x):
        for j in range(num_y):
            counts_Si[i][j] =Compute_Cell_Avarage(i,j,counts_Si)
            counts_Ox[i][j] =Compute_Cell_Avarage(i,j,counts_Ox)
            nsi = counts_Si[i][j]
            nox = counts_Ox[i][j]
            if(nsi == 0):
                ratio[i][j] = 2
            elif((nox/nsi) > 2):
                ratio[i][j] = 2
            else:
                ratio[i][j] = nox/nsi
    return counts_Si, counts_Ox, ratio


def Avarage_Hist(dat):
    tmp = np.zeros(np.shape(dat))
    nx, ny = np.shape(dat)
    for i in range(nx):
        for j in range(ny):
            dat[i][j] = Compute_Cell_Avarage(i,j,dat)

def Compute_Cell_Avarage(r, c, dat):
    idx = np.array([-1,0,1], dtype=int)
    val = 0
    for dr in idx:
        for dc in idx:
            a, b = Array_Boundary_Wrapper(r + dr, np.shape(dat)[0]), Array_Boundary_Wrapper(c + dc, np.shape(dat)[1])
            val += dat[a][b]
    return val/9


def Array_Boundary_Wrapper(x, ub):
    if(x >= 0 and x <= ub-1):
        return x
    elif(x < 0):
        return ub + x
    elif(x > ub -1):
        return x - ub
    else:
        quit("what are you doing???")



def main():
    _reg = np.array([[0,43.44], [0,27.15], [0.0,16.2]])
    _sim = rl.Read_Dump(FILE)
    dp = _plot.Dump_Plot_Map(_sim)
    dp.Plot_Map()

if __name__ =='__main__':
    test_file = "/Users/diggs/Desktop/TOPCon/out-1-9-23/SiOx-1.59.dat"
    ARGS = sys.argv
    FILE=ARGS[1]
    main()
















