import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import re,os
import matplotlib.cm as cm
import sys
from kuibit.simdir import SimDir



class Simulation:
    def __init__(self,sim_dir):
        self.proj_dir   = os.getcwd()
        self.sim_dir    = self.proj_dir+"/simulations/"+sim_dir
        self.param_file = self.proj_dir+"/parfiles/"+sim_dir+".par"
        self.metadata   = self.proj_dir+"/simulations/"+sim_dir+"/TwoPunctures.bbh"

        self.base_name = self.sim_dir.split("_")[1:]
        self.mp,  self.mm,  self.par_b, self.ex_r,\
        self.dx,  self.dy,  self.dz, \
        self.p1x, self.p1y, self.p1z,\
        self.p2x, self.p2y, self.p2z,\
        self.s1x, self.s1y, self.s1z,\
        self.s2x, self.s2y, self.s2z = self.get_params()

        self.p1 = np.sqrt(self.p1x * self.p1x +\
                          self.p1y * self.p1y +\
                          self.p1z * self.p1z  )

        self.p2 = np.sqrt(self.p2x * self.p2x +\
                          self.p2y * self.p2y +\
                          self.p2z * self.p2z  )

        self.s1 = np.sqrt(self.s1x * self.s1x +\
                          self.s1y * self.s1y +\
                          self.s1z * self.s1z  )

        self.s2 = np.sqrt(self.s2x * self.s2x +\
                          self.s2y * self.s2y +\
                          self.s2z * self.s2z  )

        self.xmesh, self.ymesh, self.umesh,\
        self.psimesh, self.bounds, self.lenx,\
        self.leny = self.load_data()

        self.yvec = self.ymesh[:,0]
        self.xvec = self.xmesh[0,:]

        self.xmin = self.bounds[0,0] ; self.xmax = self.bounds[0,1]
        self.ymin = self.bounds[1,0] ; self.ymax = self.bounds[1,1]
        self.umesh_masked   = self.mask_data("puncture_u")
        self.psimesh_masked = self.mask_data("my_psi")
        self.umax_ind,   self.umax   = self.find_max("puncture_u")
        self.psimax_ind, self.psimax = self.find_max("my_psi")


    def get_line(self,string,fl):
        file = open(fl,"r")
        for line in file:
            if re.search(string,line):
                tmp = line
                break
        file.close()
        return tmp

    def get_nvalue(self,line):
        val = line.split("=")[1]
        return float(val)

    def read_param(self,param,fl):
        return self.get_nvalue(self.get_line(param,fl))


    def get_params(self):
        mm    = self.read_param("bare-mass2",self.metadata)
        mp    = self.read_param("bare-mass1",self.metadata)
        par_b = self.read_param("position1x",self.metadata)
        ex_r  = self.read_param("excision",self.metadata)

        dx    = self.read_param("dx",self.param_file)
        dy    = self.read_param("dy",self.param_file)
        dz    = self.read_param("dz",self.param_file)

        ## Momenta
        p1x   = self.read_param("momentum1x",self.metadata)
        p1y   = self.read_param("momentum1y",self.metadata)
        p1z   = self.read_param("momentum1z",self.metadata)
        p2x   = self.read_param("momentum2x",self.metadata)
        p2y   = self.read_param("momentum2y",self.metadata)
        p2z   = self.read_param("momentum2z",self.metadata)
        ## Spins
        s1x   = self.read_param("spin1x",self.metadata)
        s1y   = self.read_param("spin1y",self.metadata)
        s1z   = self.read_param("spin1z",self.metadata)
        s2x   = self.read_param("spin2x",self.metadata)
        s2y   = self.read_param("spin2y",self.metadata)
        s2z   = self.read_param("spin2z",self.metadata)
        return mp,mm,par_b,ex_r,dx,dy,dz,p1x,p1y,p1z,p2x,p2y,p2z,s1x,s1y,s1z,s2x,s2y,s2z



    def load_data(self):
        sim = SimDir(self.sim_dir)

        u   = sim.gf.xy.fields.puncture_u[0]
        psi = sim.gf.xy.fields.my_psi[0]

        ## Use this with AMR
        #u   = u.merge_refinement_levels(resample=True)
        #psi = psi.merge_refinement_levels(resample=True)
        xmesh = u.coordinates()[0][0][0][:].T
        ymesh = u.coordinates()[1][0][0][:].T


        # filename = "/twopunctures-"+str(var)+".xy.asc-pp.txt"
        # data = pd.read_csv(self.sim_dir+filename,sep=",",comment='#',usecols=[9,10,12],names=["x","y","var"])
        # x = data['x'].values
        # y = data['y'].values

        xmin = xmesh.min()
        xmax = xmesh.max()

        ymin = ymesh.min()
        ymax = ymesh.max()

        lenx = xmesh.shape[1]
        leny = xmesh.shape[0]

        u    = u[0][0].data.T
        psi  = psi[0][0].data.T
        # var = data['var'].values

        # xmesh = x.reshape((leny,lenx))
        # ymesh = y.reshape((leny,lenx))
        # var = var.reshape((leny,lenx))
        bounds = np.array([[xmin,xmax],[ymin,ymax]])
        return xmesh,ymesh,u,psi,bounds,lenx,leny


    def mask_data(self,param):
        radius = self.ex_r
        x1d = self.xmesh.ravel()
        y1d = self.ymesh.ravel()
        if param == "puncture_u":
            var = self.umesh.copy()
        elif param == "my_psi":
            var = self.psimesh.copy()
        var1d = var.ravel()
        N   = x1d.shape[0]
        for i in range(N):
            circle = (x1d[i]+self.par_b)*(x1d[i]+self.par_b)+\
                y1d[i]*y1d[i]

            if circle < (radius*radius):
                var1d[i] = np.nan

            if x1d[i] > 0:
                var1d[i] = np.nan
        maskedmesh = var1d.reshape((self.leny,self.lenx))
        return maskedmesh

    def find_max(self,param):
        if param == "puncture_u":
            # print("puncture_u")
            var = self.umesh_masked
        elif param == "my_psi":
            # print("psi")
            var = self.psimesh_masked
        else:
            raise TypeError("Invalid choice of parameter\
            in find max")

        ind = np.unravel_index(np.nanargmax(var, axis=None), \
                               var.shape)


        return ind, var[ind]


    def psi_value_at_excision(self,hist=False):
        x1d   = self.xmesh.copy().ravel()
        y1d   = self.ymesh.copy().ravel()

        var   = self.psimesh.copy()
        var1d = var.ravel()
        N     = x1d.shape[0]
        psis  = []
        xs    = []
        ys    = []
        for i in range(N):
            circle = (x1d[i]+self.par_b)*(x1d[i]+self.par_b)+\
                y1d[i]*y1d[i]
            if circle == (self.ex_r*self.ex_r):
                psis.append(var1d[i])
                xs.append(x1d[i])
                ys.append(y1d[i])

        if self.ex_r != 0 and hist==True:
            plt.clf()
            plt.tick_params(direction="in")
            bns = int(np.sqrt(len(psis)))
            plt.hist(psis,bins=len(psis),density=True)
            plt.title(r"$\psi$ at excision sphere")
            plt.savefig("%s-hist.png" %self.sim_dir)

        return xs,ys,psis

    def plot_u_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(self.umesh, interpolation='spline36',\
                        cmap="inferno_r",origin='lower',\
                        extent=[self.xmin, self.xmax,\
                                self.ymin, self.ymax],\
                        vmin=abs(self.umesh).min(),\
                        vmax=abs(self.umesh).max())
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter([-self.par_b,self.par_b],\
                    [0,0],color="black", marker="x")
        excision_sphere = plt.Circle((-self.par_b, 0),\
                            self.ex_r,\
                            color='black',fill=False,ls="--")

        rs_radius1 = plt.Circle((-self.par_b, 0),\
                            self.mm/2,\
                            color='black',fill=False)

        rs_radius2 = plt.Circle((self.par_b, 0),\
                            self.mp/2,\
                            color='black',fill=False)

        xs, ys, psis = self.psi_value_at_excision()
        plt.scatter(xs,ys,color="green",marker="+")
        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        plt.xlim(self.xmin,self.xmax)
        plt.ylim(self.ymin,self.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r"u with Excision Radius: %.2f" %self.ex_r)
        plt.colorbar(im,orientation="horizontal")
        plt.savefig("%s-u.png" %self.sim_dir)


    def plot_psi_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(self.psimesh, interpolation='spline36',\
                        cmap="inferno_r", origin='lower',\
                        extent=[self.xmin, self.xmax,\
                                self.ymin, self.ymax],\
                        vmin=abs(self.psimesh).min(),\
                        vmax=abs(self.psimesh).max())
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter([-self.par_b,self.par_b],\
                    [0,0],color="black", marker="x")
        excision_sphere = plt.Circle((-self.par_b, 0),\
                            self.ex_r,\
                            color='black',fill=False,ls="--")

        rs_radius1 = plt.Circle((-self.par_b, 0),\
                            self.mm/2,\
                            color='black',fill=False)

        rs_radius2 = plt.Circle((self.par_b, 0),\
                            self.mp/2,\
                            color='black',fill=False)

        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        xs, ys, psis = self.psi_value_at_excision()
        plt.scatter(xs,ys,color="green",marker="+")

        plt.xlim(self.xmin,self.xmax)
        plt.ylim(self.ymin,self.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r" $\psi$ with Excision Radius: %.2f" %self.ex_r)
        plt.colorbar(im,orientation="horizontal")
        plt.savefig("%s-psi.png" %self.sim_dir)





class Error:
    def __init__(self,vnl_sim_dir):
        self.sim_dir1  = vnl_sim_dir
        self.sim_dir2  = "excision"+ vnl_sim_dir.split("vanilla")[1]
        self.vanilla   = Simulation(self.sim_dir1)
        self.excision  = Simulation(self.sim_dir2)
        self.error_u   = abs(self.vanilla.umesh -\
                             self.excision.umesh)
        self.error_u_masked = abs(self.vanilla.umesh_masked -\
                                  self.excision.umesh_masked)

        self.error_psi = abs(self.vanilla.psimesh -\
                             self.excision.psimesh)
        self.error_psi_masked = abs(self.vanilla.psimesh_masked -\
                                    self.excision.psimesh_masked)

        self.xmin      = self.vanilla.xmin
        self.xmax      = self.vanilla.xmax
        self.ymin      = self.vanilla.ymin
        self.ymax      = self.vanilla.ymax
        self.ex_r      = self.excision.ex_r
        self.indu, self.x_at_umax, self.y_at_umax, self.umax,\
            self.umaxangle = self.find_max("puncture_u")
        self.indpsi, self.x_at_psimax, self.y_at_psimax,\
            self.psimax, self.psimaxangle = self.find_max("my_psi")


    def find_max(self,param):
        if param == "puncture_u":
            # print("puncture_u")
            var = self.error_u_masked.copy()
        elif param == "my_psi":
            # print("psi")
            var = self.error_psi_masked.copy()
        else:
            raise TypeError("Invalid choice of parameter\
            in find max")

        ind = np.unravel_index(np.nanargmax(var, axis=None),\
                               var.shape)

        angle = np.arctan2(self.excision.ymesh[ind],self.excision.xmesh[ind])
        return ind, self.vanilla.xmesh[ind],self.vanilla.ymesh[ind],var[ind], angle

    def error_report(self):
        error_dict= {"x":self.x_at_psimax,\
                     "y":self.y_at_psimax,\
                     "q":self.vanilla.mp,\
                     "b":self.vanilla.par_b,\
                     "ex_r":self.ex_r,\
                     # 1st BH
                     "p1x":self.vanilla.p1x,\
                     "p1y":self.vanilla.p1y,\
                     "p1z":self.vanilla.p1z,\
                     "p1":self.vanilla.p1,\

                     "s1x":self.vanilla.s1x,\
                     "s1y":self.vanilla.s1y,\
                     "s1z":self.vanilla.s1z,\
                     "s1":self.vanilla.s1,\



                     # 2nd BH
                     "p2x":self.vanilla.p2x,\
                     "p2y":self.vanilla.p2y,\
                     "p2z":self.vanilla.p2z,\
                     "p2":self.vanilla.p2,\

                     "s2x":self.vanilla.s2x,\
                     "s2y":self.vanilla.s2y,\
                     "s2z":self.vanilla.s2z,\
                     "s2":self.vanilla.s2,\


                     "angle":self.psimaxangle,\
                     "max_error_psi":self.psimax}
        return pd.DataFrame.from_dict([error_dict])

    def plot_error_u_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(self.error_u, interpolation='spline36',\
                        cmap="inferno_r", origin='lower',\
                        extent=[self.xmin, self.xmax,\
                                self.ymin, self.ymax],\
                        vmin=abs(self.error_u).min(),\
                        vmax=abs(self.error_u).max())
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter([-self.vanilla.par_b,self.vanilla.par_b],\
                    [0,0],color="black", marker="x")

        excision_sphere = plt.Circle((-self.excision.par_b, 0),\
                            self.excision.ex_r,\
                            color='black',fill=False,ls="--")

        rs_radius1 = plt.Circle((-self.excision.par_b, 0),\
                            self.vanilla.mm/2,\
                            color='black',fill=False)

        rs_radius2 = plt.Circle((self.excision.par_b, 0),\
                            self.vanilla.mp/2,\
                            color='black',fill=False)

        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        plt.scatter([self.x_at_umax],[self.y_at_umax],color="green",marker="+")



        plt.xlim(self.vanilla.xmin,self.vanilla.xmax)
        plt.ylim(self.vanilla.ymin,self.vanilla.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r"Error_u with Excision Radius: %.2f" %self.ex_r)
        plt.colorbar(im,orientation="horizontal")
        plt.savefig("%s-error-u.png" %self.vanilla.sim_dir)
        plt.savefig("%s-error-u.png" %self.excision.sim_dir)

    def plot_error_psi_xy(self):
        plt.clf()
        plt.tick_params(direction="in")
        im = plt.imshow(self.error_psi, interpolation='spline36',\
                        cmap="inferno_r", origin='lower',\
                        extent=[self.xmin, self.xmax,\
                                self.ymin, self.ymax],\
                        vmin=abs(self.error_psi).min(),\
                        vmax=abs(self.error_psi).max())
        # norm= colors.Normalize(vmin=u.min(), vmax=u.max())
        plt.scatter([-self.vanilla.par_b,self.vanilla.par_b],\
                    [0,0],color="black", marker="x")
        excision_sphere = plt.Circle((-self.excision.par_b, 0),\
                            self.excision.ex_r,\
                            color='black',fill=False,ls="--")

        rs_radius1 = plt.Circle((-self.excision.par_b, 0),\
                            self.vanilla.mm/2,\
                            color='black',fill=False)

        rs_radius2 = plt.Circle((self.excision.par_b, 0),\
                            self.vanilla.mp/2,\
                            color='black',fill=False)

        plt.gcf().gca().add_artist(excision_sphere)
        plt.gcf().gca().add_artist(rs_radius1)
        plt.gcf().gca().add_artist(rs_radius2)

        plt.scatter([self.x_at_psimax],[self.y_at_psimax],color="red",marker="+")


        plt.xlim(self.vanilla.xmin,self.vanilla.xmax)
        plt.ylim(self.vanilla.ymin,self.vanilla.ymax)
        plt.xlabel("x/M")
        plt.ylabel("y/M")
        plt.title(r"Error_$\psi$ with Excision Radius: %.2f" %self.ex_r)
        plt.colorbar(im,orientation="horizontal")
        plt.savefig("%s-error-psi.png" %self.vanilla.sim_dir)
        plt.savefig("%s-error-psi.png" %self.excision.sim_dir)








# vnl = Simulation("vanilla")
# exc = Simulation("excision")
# x,y,psi = exc.psi_value_at_excision()
# print(np.mean(psi))
# print(exc.mm/(2*exc.ex_r))

# vnl.plot_psi_xy()
# vnl.plot_u_xy()

# exc.plot_psi_xy()
# exc.plot_u_xy()

# error = Error()
# error.plot_error_u_xy()
# error.plot_error_psi_xy()

class BlackHole:
    def __init__(self,sim_dir):
        proj_dir        = "/work/stamatis.vretinaris/imri-id/"
        self.sim_dir    = proj_dir+"simulations/dummy_"+sim_dir+"/"+sim_dir
        self.param_file = proj_dir+"parfiles/"+sim_dir+".par"
        self.metadata   = proj_dir+"simulations/dummy_"+sim_dir+"/"+sim_dir+"/TwoPunctures.bbh"

        self.spinx     = 0.0
        self.spiny     = 0.0
        self.spinz     = 0.0
        self.momentumx = 0.0
        self.momentumy = 0.0
        self.momentumz = 0.0
        self.mass      = 1.0

    def get_line(self,string,fl):
        file = open(fl,"r")
        for line in file:
            if re.search(string,line):
                tmp = line
                break
        file.close()
        return tmp

    def get_nvalue(self,line):
        val = line.split("=")[1]
        return float(val)

    def read_param(self,param,fl):
        return self.get_nvalue(self.get_line(param,fl))


    def get_params(self):
        mm    = self.read_param("bare-mass2",self.metadata)
        mp    = self.read_param("bare-mass1",self.metadata)
        par_b = self.read_param("position1x",self.metadata)
        ex_r  = self.read_param("excision",self.metadata)

        dx    = self.read_param("dx",self.param_file)
        dy    = self.read_param("dy",self.param_file)
        dz    = self.read_param("dz",self.param_file)
        ## Momenta
        p1x   = self.read_param("momentum1x",self.metadata)
        p1y   = self.read_param("momentum1y",self.metadata)
        p1z   = self.read_param("momentum1z",self.metadata)
        p2x   = self.read_param("momentum2x",self.metadata)
        p2y   = self.read_param("momentum2y",self.metadata)
        p2z   = self.read_param("momentum2z",self.metadata)
        ## Spins
        s1x   = self.read_param("spin1x",self.metadata)
        s1y   = self.read_param("spin1y",self.metadata)
        s1z   = self.read_param("spin1z",self.metadata)
        s2x   = self.read_param("spin2z",self.metadata)
        s2y   = self.read_param("spin2y",self.metadata)
        s2z   = self.read_param("spin2z",self.metadata)

        return mp,mm,par_b,py,ex_r,dx,dy,dz
