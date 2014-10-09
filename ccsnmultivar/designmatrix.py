# design matrix class and methods
import re
import warnings
import numpy as np
import csv

class DesignMatrix(object):
    def __init__(self,path_to_parameterfile,formula):
        self._path = path_to_parameterfile
        self.formula = formula
        # call functions to fit on instantiation
        self._load_paramfile()
        self.fit_transform()

        # print warning about degrees of freedom
        if np.shape(self._X)[0] <= np.shape(self._X)[1]:
            # if n < p
            warnings.warn("Number of columns > number of waveforms.  Hotellings T2 not defined")

    def _load_paramfile(self):
        """
        Loads up the parameter file using supplied directory from self._path

        Parameter files must be:
        - in .csv or .dat format
        - first row is the header
        - first column is a string giving the waveform name
        - other columns are values of physical parameter for each waveform
        """
        # open from parameter path
        param_list = list(csv.reader(open(self._path,"rb")))
        # delete empty elements (if any)
        self.param_file = [x for x in param_list if x != []]

    def fit_transform(self):
        """
        Transform raw parameters to design matrix.
        Returns: X_dict, wave_names
            X_dict[column] = np.array(values)
            wavenames in correct order to correspond to rows in np.array(values)
        """
        param_file = self.param_file
        formula = self.formula
        # make dict of [wavenames] = raw_params
        name_list = []
        param_list = []
        # get header names for each param (names of param_list columns)
        param_colnames = param_file[0][1:] # 0th element is "Name" or "Wavename"
        # start from 1.  (row 0 is the header)
        for i in np.arange(1, len(param_file)):
            name_list.append(param_file[i][0])
            param_list.append(param_file[i][1:])
        # remove ' ' blank spaces from param_list
        param_list = [[x.strip() for x in y] for y in param_list]
        param_dict = {}
        # double for loop. i loops through param_colnames, j loops thru param values per wave
        for i in np.arange(0, len(param_colnames)):
            param_dict[param_colnames[i]] = []
            for j in np.arange(0,len(name_list)):
                param_dict[param_colnames[i]].append(param_list[j][i])
        # now we have param_dict, and name_list
        # parse formula
        formula_dict, inter_list = _parse_formula(formula)
        # turn instructions and raw parameters into design matrix
        X_final,col_names = _encode_design_matrix(formula_dict,inter_list,param_dict)
        # add intercept term to X_final, add 'Intercept' to first in col_names
        warnings.warn("Design Matrix includes Intercept Term")
        col_names.insert(0,'Intercept')
        X_final = np.concatenate([np.ones((len(name_list),1)),X_final],axis=1)

        # put design matrix, row names, and column names, into class variables
        self._X = X_final
        self._X_col_names = col_names
        self._row_names   = name_list

    def get_columnnames(self):
        return self._X_col_names
    def get_rownames(self):
        return self._row_names
    def get_X(self):
        return np.asarray(self._X)
    def get_formula(self):
        return self._formula


def _parse_formula(formula):
    """
    Parse formula into a dictionary
    formula_dict[variable_name] = [encoding, dropped_name]
    Parse interactions into a list
    inter_list = [[A,B], [A,C], [A,B,C]]

    formula = "A + beta + A*beta | Dev(A,drop=1), Poly(beta,degree=3)"
    """
    #TODO: DEAL WITH BSPLINE HAS MULTIPLE ARGUEMENTS (NOT JUST ONE)

    # break formula apart from encoding instructions
    formula,instr = formula.replace(' ','').split('|')
    # break formula apart at + sign
    formula_terms = formula.split('+') 
    # examine the instructions term, first split by ),
    instr = instr.split('),')
    # elements in the instr list will match the number of non interaction 
    #      elements in formula
    # go through each formula term, making a dictionary whose key is variable name
    formula_dict = {}
    encoding     = []
    other_arg    = []
    inter_list   = []
    for term in formula_terms:
        if "*" in term:
            # then this is an interaction term, make 'list of lists'
            term = term.split('*')
            inter_list.append(term)
        else:
            # then this is not an interaction term, make blank dictionary
            formula_dict[term] = [encoding, other_arg]
    # loop through instructions, parse each term
    for term in instr:
        # remove punctuation in term
        term = re.sub('[()]','',term)
        # check for each encoding type
        if "Dev" in term:
            # remove first three letters (Dev)
            term = term[3:]
            # split on comma
            var_name,arg = term.split(',')
            # split on equals sign in arg, take part after
            drop_name = arg.split('=')[1]
            # put var_name and drop_name into proper key in formula_dict
            formula_dict[var_name] = ["Dev",drop_name]
        elif "Dum" in term:
            # remove first three letters (Dum)
            term = term[3:]
            # split on comma
            var_name,arg = term.split(',')
            # split on equals sign in arg, take part after
            ref_name = arg.split('=')[1]
            # put var_name and drop_name into proper key in formula_dict
            formula_dict[var_name] = ["Dum",ref_name]
        elif "Poly" in term:
            # remove first four letters (Poly)
            term = term[4:]
            # split on comma
            var_name,arg = term.split(',')
            # split on equals sign in arg, take part after
            degree = arg.split('=')[1]
            # put var_name and drop_name into proper key in formula_dict
            formula_dict[var_name] = ["Poly",degree]
        else:
            raise Exception("Unknown Encoding")
    return formula_dict,inter_list



def _encode_design_matrix(formula_dict, inter_list, param_dict):
    """
    1. Deviation encoding
    2. Dummy encoding
    4. Chebyshev Polynomial encoding
    6. Interaction between any of the above types
    """
    X = []
    col_names = []
    X_dict = {}
    Xcol_dict = {}
    # first, replace param_dict[key] = values, with param_dict[key] = dmatrix
    for key in formula_dict:
        encoding,arg = formula_dict[key]
        if 'Dev' in encoding:
            # make deviation encoded design matrix
            drop_name = arg
            # encode
            X_sub,colnames_sub = _dev_encode(param_dict,drop_name,key)
            # additionally, store in dictionary for use by interactions
            X_dict[key] = X_sub
            Xcol_dict[key] = colnames_sub
        elif 'Dum' in encoding:
            # make dummy variable encoding design mat
            ref_name = arg
            X_sub,colnames_sub = _dum_encode(param_dict,ref_name,key)
            # additionally, store in dictionary for use by interactions
            X_dict[key] = X_sub
            Xcol_dict[key] = colnames_sub
        elif 'Poly' in encoding:
            # make polynomial encoding design mat
            degree = arg
            X_sub,colnames_sub = _poly_encode(param_dict,degree,key)
            # additionally, store in dictionary for use by interactions
            X_dict[key] = X_sub
            Xcol_dict[key] = colnames_sub
        else:
            raise Exception("Encoding name error")
        # update with each new encoding
        X.append(X_sub)
        col_names.append(colnames_sub)
    # now compute interaction designmatrices
    for interaction in inter_list:
        if len(interaction) >= 3:
            raise Exception("Doesn't allow 4-way or higher interaction terms")
        elif len(interaction) == 3:
            # there are three interaction terms (A*B*C)
            # pull them out into X_i and names_i matrices and list of strings
            X1 = X_dict[interaction[0]]
            names_1 = Xcol_dict[interaction[0]]
            X2 = X_dict[interaction[1]]
            names_2 = Xcol_dict[interaction[1]]
            X3 = X_dict[interaction[2]]
            names_3 = Xcol_dict[interaction[2]]

            X_int = []
            names_int = []
            for i in np.arange(0,X1.shape[1]):
                for j in np.arange(0,X2.shape[1]):
                    for k in np.arange(0,X3.shape[1]):
                        X_int.append(X1[:,i]*X2[:,j]*X3[:,k])
                        col_names.append(names_1[i] + "*" + names_2[j] + "*" + names_3[k])
            # make X_int from lists to np array
            X_int = np.array(X_int).T

        elif len(interaction) == 2:
            # there are two interaction terms (A*B)
            # pull them out into X_i and names_i matrices and list of strings
            X1 = X_dict[interaction[0]]
            names_1 = Xcol_dict[interaction[0]]
            X2 = X_dict[interaction[1]]
            names_2 = Xcol_dict[interaction[1]]

            X_int = []
            names_int = []
            for i in np.arange(0,X1.shape[1]):
                for j in np.arange(0,X2.shape[1]):
                    X_int.append(X1[:,i]*X2[:,j])
                    col_names.append(names_1[i] + "*" + names_2[j])
            X_int = np.array(X_int).T
        else:
            raise Exception("Error while evaluating meaning of interaction term")
        # now concatenate X_int and names_int for the next interaction term (if any)
        X.append(X_int)
        col_names.append(names_int)
        # convert X from list of design matrices to design matrix
        X = np.concatenate(X,axis=1)
    return X, col_names


def _dev_encode(param_dict,drop_name,param_name):
    # param_dict[col_name] = np.array(wave_values)
    values = param_dict[param_name]
    # treat as string, map each unique value to an integer
    unique_values = []
    [unique_values.append(i) for i in values if not unique_values.count(i)]
    # warning if number of unique values is large compared to number of waveforms
    if len(unique_values) > 0.5*len(values):
        warning.warn("Many unique values in parameter %s for a Deviation encoding" & param_name)
    # make deviation encoding instruction matrix
    I = np.eye(len(unique_values))
    # remove column of I for drop_name
    drop_idx = unique_values.index(drop_name)
    I = np.delete(I,drop_idx,1)
    # make row of I for drop_name index
    I[drop_idx,:] = -1
    # make names out of values
    col_names = []
    for i in np.arange(0,len(unique_values)):
        col_names.append(param_name+":["+unique_values[i]+" - "+"mu"+"]")
    # remove column name for drop_name
    col_names.remove(param_name + ":[" + drop_name + " - " + "mu" + "]")
    # use instruction matrix I to map parameter values to design matrix
    X = np.empty((len(values),len(col_names)))
    for i in np.arange(0,len(values)):
        idx = unique_values.index(values[i])
        X[i,:] = I[idx,:]
    return X, col_names


def _dum_encode(param_dict,drop_name,param_name):
    # param_dict[col_name] = np.array(wave_values)
    values = param_dict[param_name]
    # treat as string, map each unique value to an integer
    unique_values = []
    [unique_values.append(i) for i in values if not unique_values.count(i)]
    # warning if number of unique values is large compared to number of waveforms
    if len(unique_values) > 0.5*len(values):
        warning.warn("Many unique values in parameter %s for a Deviation encoding" & param_name)
    # make deviation encoding instruction matrix
    I = np.eye(len(unique_values))
    # remove column of I for drop_name
    drop_idx = unique_values.index(drop_name)
    I = np.delete(I,drop_idx,1)
    # make names out of values
    col_names = []
    for i in np.arange(0,len(unique_values)):
        col_names.append(param_name+":["+unique_values[i]+" - " + drop_name + "]")
    # remove column name for drop_name
    col_names.remove(param_name + ":[" + drop_name + " - " + drop_name + "]")
    # use instruction matrix I to map parameter values to design matrix
    X = np.empty((len(values),len(col_names)))
    for i in np.arange(0,len(values)):
        idx = unique_values.index(values[i])
        X[i,:] = I[idx,:]
    return X, col_names



def _poly_encode(param_dict,degree,param_name):
    x = np.array(map(float,param_dict[param_name]))
    degree = int(degree)
    if len(np.unique(x)) < 5:
        warning.warn("Not many unique values of parameter %s for Poly encoding" & param_name)
    # generate chebyshev polynomials
    n = degree
    m = len(x)
    # generate the z variable as a mapping of your x data into [-1,1]
    z = ((x - .99*min(x))-(.99*max(x)-x))/(.99*max(x) - .99*min(x))
    X = np.empty([m,n+1])
    X[:,0] = np.ones([m,])
    if n >= 1:
        X[:,1] = z;
    if n >= 2:
        for k in np.arange(2, n+1):
            X[:,k] = 2.*z*X[:,k-1] - X[:,k-2]
    # remove first column of A (never keep intercept)
    X = X[:,1:]
    # generate names for each column
    col_names = []
    for i in np.arange(0,degree):
        col_names.append(param_name + "^" + str(i+1))
    return X, col_names

def _bs_encode(param_dict,bsplineparameters, param_name):
    print "not implemented yet"
