import scipy.sparse as sp
import numpy as np
import matplotlib.pylab as plt
import os, sys
import platform

#if platform.system() == 'Windows':
#if cpp_solver=="EMPY":
    #sys.path.append(r'C:\EMSolution\EMSolPy3\x64\Release') 
sys.path.append(r'..\bin\Release') 
import EMPY_Solver
#else:
sys.path.append(r'..\bin\Release')
import SparseSolvPy
    
class MatrixSolver:
    def __init__(self, matrix):
        self.mat=matrix
    
    def AddCoupling(self,cvecs, amat, ccomplex):
        matrix=self.mat
        dim=matrix.shape[0]
        #print("dim=", dim)
        rows, cols = matrix.nonzero()
        vals = matrix[rows, cols]
        vals = np.ravel(vals)
        size=vals.size 
        #size=matrix.nnz
    
        nadd=len(cvecs)
        non_zeros=0
        for n in range(nadd):
            m = np.count_nonzero(cvecs[n])
            if m> non_zeros: non_zeros=m
          
        print("size of matrix= ", size)
        sizep=size +non_zeros*2*nadd +nadd*nadd
        new_rows= np.zeros(sizep, dtype=int)
        new_cols= np.zeros(sizep, dtype=int)
        if ccomplex==False:
            new_vals= np.zeros(sizep)
        else:
             new_vals= np.zeros(sizep, dtype=complex)
            
        for i in range(size):
            new_rows[i] =rows[i]
            new_cols[i] =cols[i]
            new_vals[i] =vals[i]
        rows=None
        cols=None
    
        k=size
        vals=None
    
        for n in range(nadd):
            fcut = cvecs[n]
            r=dim+n
            for col in range(dim):
                v=fcut[col]
                if v !=0:
                    new_rows[k]=col
                    new_cols[k]=r
                    new_vals[k]=v
                    k=k+1
                    new_rows[k]=r
                    new_cols[k]=col
                    new_vals[k]=v
                    k=k+1
                
        for n in range(nadd):             
            for n2 in range(nadd):        
                new_rows[k]=dim+n
                new_cols[k]=dim+n2
                new_vals[k]=amat[n][n2]
                k=k+1
 
        new_a=sp.csr_matrix((new_vals, (new_rows, new_cols)), shape=(dim+nadd, dim+nadd))
        new_rows=None
        new_cols=None
        new_vals=None
        rows, cols = new_a.nonzero()
        vals = new_a[rows, cols]
        vals = np.ravel(vals)
    
        #print("rows=", rows)
        #print("cols=", cols)
        #print("vals=", vals)
        
        return new_a;
    
    @classmethod
    def iccg_solve(cls,fes, gf, A, Bvec,  **kwargs):
        import time
        print("enter iccg_solve")
        default_values = {"tol": 1.E-10, "accel_factor":1.1, "max_iter":1000, "complex":False, 
                          "divfac":10., "diviter":10,
                          "logplot":False, "cpp_solver":"EMPY",
                          "scaling":True}
        default_values.update(kwargs)
        tol=default_values["tol"]
        accel_factor=default_values["accel_factor"]
        max_iter=default_values["max_iter"]
        ccomplex=default_values["complex"]
        divfac=default_values["divfac"]
        diviter=default_values["diviter"]
        logplot=default_values["logplot"]
        cpp_solver=default_values["cpp_solver"]
        scaling=default_values["scaling"]
    
        asci = sp.csr_matrix (A.mat.CSR())
        #    A=None
        Acut = asci[:,fes.FreeDofs()][fes.FreeDofs(),:]
        asci=None
        fcut = np.array(Bvec)[fes.FreeDofs()]
        #print("in iccg_solve")
        #print("fcut=", fcut)
        ucut = np.array(Bvec, copy=True)[fes.FreeDofs()]
        rows, cols = Acut.nonzero()
        vals = Acut[rows, cols]
        #Acut=None
        vals = np.ravel(vals)
        dim=fcut.size
        size= (len(rows)-dim)/2-dim
        print('Dof=',dim, '   Nonzeros=', len(rows))
        if dim ==0: return gf

        #if platform.system() == 'Windows':
        if cpp_solver=="EMPY":
            if ccomplex==True:
                solver=EMPY_Solver.EMPY_CSolver();
            else:
                solver=EMPY_Solver.EMPY_Solver();

            #solver.Write(dim, len(rows), rows, cols, vals, fcut)
            solver.SetMatrix(dim, len(rows), rows, cols, vals)
            rows=None
            cols=None
            vals=None
            solver.SetScaling(scaling);
            #solver.SetMethod("ICCG")
            solver.SetEps(tol);
            solver.SetShiftParameter(accel_factor);
            solver.SetDivCriterion(divfac,  diviter)
            
            start_time = time.perf_counter()           
            ucut=solver.Solve(fcut, ucut)
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            
            shift=solver.GetShiftParameter()
            print("shift parameter=", shift)
            log1 = solver.GetResidualLog()
            #print(log1)
            its=solver.GetMinimumResidual()
            if len(log1) !=0: print("minimum residual=", min(log1), " at iteraions: ", its)

            #print("log1  ",log1)
        else:
            if ccomplex==True: 
                mat = SparseSolvPy.SparseMatC(len(rows), rows, cols, vals)
            else:
                mat = SparseSolvPy.SparseMat(len(rows), rows, cols, vals) 
            """  
            import pickle
            my_list=[dim, len(rows), rows, cols, vals, fcut]
            with open('my_list.pkl', 'wb') as f:
                pickle.dump(my_list, f)

            with open('my_list.pkl', 'rb') as f:
                my_list_read = pickle.load(f)
   
            dim=my_list_read[0]
            lenrows=my_list_read[1]
            rows=my_list_read[2]
            cols=my_list_read[3]
            vals=my_list_read[4]
            fcut=my_list_read[5]
            print("dim=", dim, "lenrows=", lenrows, "rows=", rows);
            if ccomplex==True: 
                mat = SparseSolvPy.SparseMatC(len(rows), rows, cols, vals)
            else:
                mat = SparseSolvPy.SparseMat(len(rows), rows, cols, vals) 
            """
            rows=None
            cols=None
            vals=None
            if accel_factor <= 0: accel_factor=1.1
            solver = SparseSolvPy.MatSolvers()
            solver.setSaveBest(True)
            solver.setSaveLog(True)
            solver.setDiagScale(False)
            solver.setDirvegeType(1)
            solver.setBadDivCount(10)
            
            start_time = time.perf_counter()           
            ucut=solver.Solve(fcut, ucut)
            solver.solveICCG_py(len(fcut), tol, max_iter, accel_factor, mat, fcut, ucut, True)
            elapsed_time = end_time - start_time
            log1 = solver.getResidualLog_py()

        #print("log1  ",log1)
        if logplot==True:
            plt.plot(range(len(log1)), log1)    
            plt.yscale('log')
            plt.show(block=False)

        np.array(gf.vec.FV(), copy=False)[fes.FreeDofs()] += ucut
    
        result = Acut.dot(ucut) - fcut 
        #norm = np.linalg.norm(result)/np.linalg.norm(fcut)
        print("ICCG calculation time (sec):", elapsed_time)
        #power=np.dot(fcut, ucut)
        #print("power= ", power)
    
        ucut=None
        mat=None
        fcut=None

        #log1min.append(min(log1))
        return gf #, power
    
    @classmethod
    def SolveCoupled(cls, fes, A, c, a, f, g, **kwargs):
        default_values = {"tol": 1.E-10, "accel_factor":1.1, "max_iter":1000}
        default_values.update(kwargs)
        tol=default_values["tol"]
        accel_factor=default_values["accel_factor"]
        max_iter=default_values["max_iter"]
        voltage=g
        asci = sp.csr_matrix (A.mat.CSR())
        print("shape=", asci.shape)
        Acut = asci[:,fes.FreeDofs()][fes.FreeDofs(),:]
        asci=None
        dim=Acut.shape[0]
        nadd=len(c)
        print("dim=", dim, "  nadd=", nadd)

        ms=MatrixSolver(matrix)
        new_a=m.AddCoupling(c, a)

        
        rows, cols = new_a.nonzero()
        vals = new_a[rows, cols]
        vals = np.ravel(vals)
    
        mat = SparseSolvPy.SparseMat(len(rows), rows, cols, vals)    
        solver = SparseSolvPy.MatSolvers()
        solver.setSaveBest(True)
        solver.setSaveLog(True)
        solver.setDiagScale(False)
        solver.setDirvegeType(1)
        solver.setBadDivCount(10)
        solver.setBadDivVal(10.0)
        fcut = np.zeros(dim)
        loop="inedge"
        if f== None:
            f=np.zeros(dim)
        else:    
            f= np.array(f)[fes.FreeDofs()]
        f=np.append(f, voltage)
        u = np.zeros(dim+nadd)
        solver.solveICCG_py(dim+1, tol, max_iter, accel_factor, mat, f, u, True)
        x =u[:-1]
        y=u[dim]
        u=None

        mat=None
        log1 = solver.getResidualLog_py()
    #    print(log1)
        plt.plot(range(len(log1)), log1)    
        plt.yscale('log')
        plt.show(block=False)  
        print("min:", min(log1))
    
        return x, y

    @classmethod
    def SolveCoupled2(cls, fes, A, c, a, f, g, **kwargs):            
        default_values = {"tol": 1.E-10, "accel_factor":1.1, "max_iter":1000, "complex":False, "logplot":False}
        default_values.update(kwargs)
        tol=default_values["tol"]
        accel_factor=default_values["accel_factor"]
        max_iter=default_values["max_iter"]
        ccomplex=default_values["complex"]
        logplot=default_values["logplot"]

        #nloop=len(c)
        asci = sp.csr_matrix (A.mat.CSR())
        #print("shape=", asci.shape)
        #print("in iccg_solve0")
        Acut = asci[:,fes.FreeDofs()][fes.FreeDofs(),:]
        #print("in iccg_solve")
        #print("Acut=", Acut)
        asci=None
        dim=Acut.shape[0]
        #print("dim=", dim)
        nadd=len(c)
        cvecs=[]
        for n in range(nadd):
            fcut = np.array(c[n].vec.FV())[fes.FreeDofs()]
            cvecs.append(fcut)
        print("dim=", dim, "  nadd=", nadd)
        ms=MatrixSolver(Acut)
        new_a=ms.AddCoupling(cvecs, a, ccomplex)

        rows, cols = new_a.nonzero()
        vals = new_a[rows, cols]
        vals = np.ravel(vals)
    
        #print("rows=", rows)
        #print("cols=", cols)
        #print("vals=", vals)

        fcut = np.zeros(dim)
        if f== None:
            f=np.zeros(dim)
        else:    
            f= np.array(f)[fes.FreeDofs()]
        voltage=g
        for n in range(nadd):
            val=0
            if voltage!=None: val=voltage[n]
            f=np.append(f, val)
        u = np.zeros(dim+nadd)
        
        if platform.system() == 'Windows': 
            if ccomplex==True:
                solver=EMPY_Solver.EMPY_CSolver();
            else:
                solver=EMPY_Solver.EMPY_Solver();
           
            solver.SetMatrix(dim+nadd, len(rows), rows, cols, vals)
   
            rows=None
            cols=None
            vals=None
            #solver.SetScaling(True);
            solver.SetEps(tol);
            solver.SetShiftParameter(accel_factor);
            solver.SetDivCriterion(10.,  10)
 
            u=solver.Solve(f, u);
            shift=solver.GetShiftParameter()
            print("shift parameter=", shift)
            log1 = solver.GetResidualLog()
            its=solver.GetMinimumResidual()
            print("minimum residual=", min(log1), " at iteraions: ", its)
            
            if logplot==True:
                plt.plot(range(len(log1)), log1)    
                plt.yscale('log')
                plt.show(block=False)

        else:
            if ccomplex==True: 
                mat = SparseSolvPy.SparseMatC(len(rows), rows, cols, vals)
            else:
                mat = SparseSolvPy.SparseMat(len(rows), rows, cols, vals) 
            rows=None
            cols=None
            vals=None  
            solver = SparseSolvPy.MatSolvers()
            solver.setSaveBest(True)
            solver.setSaveLog(True)
            solver.setDiagScale(False)
            solver.setDirvegeType(1)
            solver.setBadDivCount(10)
            solver.setBadDivVal(10.0)
            solver.solveICCG_py(dim+nadd, tol, max_iter, accel_factor, mat, f, u, True)
            
            log1 = solver.getResidualLog_py()
            #print(log1)
            plt.plot(range(len(log1)), log1)    
            plt.yscale('log')
            plt.show(block=False)  
            print("min:", min(log1))

        x=u
        if nadd !=0: x =u[:-nadd]
        y=[]
        for n in range(nadd):
            y.append(u[dim+n])
        #u=None
        mat=None
        #print("new_a=", new_a)
        #print("u=", u)
        if ccomplex==False:
            result = new_a.dot(np.array(u, dtype=float)) - f 
        else:
            result = new_a.dot(np.array(u, dtype=complex)) - f 
        #print("np.linalg.norm(result)=", np.linalg.norm(result))
        #print("np.linalg.norm(f)=",np.linalg.norm(f))
        norm = np.linalg.norm(result)/np.linalg.norm(f)
        print("結果のノルム:", norm)  
     
        return x, y

