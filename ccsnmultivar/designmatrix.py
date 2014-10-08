# design matrix class and methods
import patsy as pt
import pandas as pd
import re
import tabulate as tabulate
import warnings
import numpy as np

class DesignMatrix(object):
    def __init__(self,path_to_parameterfile,formula):
        self._formula        = formula
        self._parameter_df   = paramfile_to_df(path_to_parameterfile)
        self._X_df           = formula_to_design_matrix(self._formula,self._parameter_df)
        if len(self._X_df) < len(self._X_df.columns.tolist()):
            # if n < p
            warnings.warn("Number of columns > number of waveforms.  Hotellings T2 not defined")
            # dof is degrees of freedom
            self._X_dof = np.nan
        else:
            self._X_dof =  len(self._X_df) - len(self._X_df.columns.tolist())

    def get_dof(self):
        return self._X_dof
    def get_columnnames(self):
        return self._X_df.columns.tolist()
    def get_X(self):
        return np.asarray(self._X_df)
    def get_formula(self):
        return self._formula

def paramfile_to_df(path_to_parameterfile):
    """
    Parameter files must be:
    - in .csv or .dat format
    - first row is the header
    - first column is a string giving the waveform name
    - other columns are values of physical parameter for each waveform
    """
    # TODO make sure there are no spaces in header names (or anywhere)
    param_df = pd.read_csv(path_to_parameterfile)
    param_df.dropna(how="all",inplace=True)
    return param_df

def formula_to_design_matrix(formula,data_frame):
    formula_dict,inter_list = parse_formula(formula)
    Xdf = encode_design_matrix(formula_dict,inter_list,data_frame)
    return Xdf


def parse_formula(formula):
    """
    Parse formula into a dictionary
    formula_dict[variable_name] = [encoding, dropped_name]
    Parse interactions into a list
    inter_list = [A:B, A:C, A:B:C]

    formula = "A + beta + A*beta | Dev(A,drop=1), Poly(beta,degree=3)"
    """
    #TODO: DEAL WITH BSPLINE HAS MULTIPLE ARGUEMENTS (NOT JUST ONE)
    #TODO: error checking (interaction A*A),
    #TODO: deal with situation when no interaction terms are passed

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
    print formula_dict, inter_list
    return formula_dict,inter_list



def encode_design_matrix(formula_dict, inter_list, data_frame):
    """
    1. Deviation encoding
    2. Dummy encoding
    4. Chebyshev Polynomial encoding    5. Spline
    6. Interaction between any of the above two types
    """
    # TODO errors when there is no interaction in the formula

    # first make intercept term
    print formula_dict.keys()
    Xdf_dict = {}
    for key in formula_dict:
        encoding,arg = formula_dict[key]
        if 'Dev' in encoding:
            # make deviation encoded design matrix
            drop_name = arg
            X,column_names = dev_encode(data_frame,drop_name,key)
            # convert X and column names to data frame
            Xdf_dict[key] = pd.DataFrame(X,columns=column_names)
        elif 'Dum' in encoding:
            # make dummy variable encoding design mat
            ref_name = arg
            X,column_names = dum_encode(data_frame,ref_name,key)
            # convert X and column names to data frame
            Xdf_dict[key] = pd.DataFrame(X,columns=column_names)
        elif 'Poly' in encoding:
            # make polynomial encoding design mat
            degree = arg
            X,column_names = poly_encode(data_frame,degree,key)
            # convert X and column names to data frame
            Xdf_dict[key] = pd.DataFrame(X,columns=column_names)
        else:
            raise Exception("Encoding name error")

    # now compute interaction dataframes
    list_of_new_inter_dfs = []
    for interaction in inter_list:
        Xdf_list = []
        for variable in interaction:
            Xdf_list.append(Xdf_dict[variable])
        # now Xdf_list is a list where each element is a df of the design 
        #    matrix involved in the interaction
        # TODO: THIS ONLY WORKS WITH TWO-WAY INTERACTIONS!
        warnings.warn("Warning: Only able to do two-way interactions!")
        new_inter_dict = {}
        for i in np.arange(0,len(Xdf_list[0].columns)):
            for j in np.arange(0,len(Xdf_list[1].columns)):
                inter_name_i = Xdf_list[0].columns.tolist()[i]
                inter_name_j = Xdf_list[1].columns.tolist()[j]
                inter_nameij = inter_name_i + "*" + inter_name_j
                column_i     = np.array(Xdf_list[0].iloc[:,i])
                column_j     = np.array(Xdf_list[1].iloc[:,j])
                new_column   = column_i*column_j
                new_inter_dict[inter_nameij] = new_column
        list_of_new_inter_dfs.append(pd.DataFrame(new_inter_dict))

    # concatenate together the list of df dictionaries, and interactiondf dictionaries
    Xdf_int = pd.concat(list_of_new_inter_dfs,axis=1)
    Xdf     = pd.concat(Xdf_list, axis=1)
    Xdf     = pd.concat([Xdf,Xdf_int],axis=1)
    # add intercept to df
    X_intercept = pd.DataFrame(np.ones(len(Xdf)), columns = ['Intercept'])
    Xdf     = pd.concat([X_intercept, Xdf],axis=1)
    return Xdf




def dev_encode(data_frame, drop_name, variable_name):
    patsy_formula = "C(" + variable_name + ", Sum(omit = " + str(drop_name) + "))"
    ptmat = pt.dmatrix(patsy_formula,data_frame)
    X = np.asarray(ptmat)
    # remove intercept
    X = X[:,1:]

    # make column names
    patsy_names = ptmat.design_info.column_names[1:]
    col_names = []
    for name in patsy_names:
        # in this encoding setting, the variable name is 4th char from end
        var_name = name[-4]
        col_names.append("["+variable_name + var_name + " - " + "mu" + "]")
    return X, col_names

def dum_encode(data_frame, ref_name, variable_name):
    patsy_formula = "C(" + variable_name + ", Treatment(reference = " + str(ref_name) + "))"
    ptmat = pt.dmatrix(patsy_formula,data_frame)
    X = np.asarray(ptmat)
    # remove intercept
    X = X[:,1:]

    # make column names
    patsy_names = ptmat.design_info.column_names[1:]
    col_names = []
    var_type = variable_name
    ref_val  = ref_name
    for name in patsy_names:
        # in this encoding setting, the comparison name is 4th char from end
        var_val = name[-4]
        col_names.append("["+var_type+var_val + " - " + var_type+ref_val+"]")
    return X, col_names






def poly_encode(data_frame,degree,variable_name):
    x = np.array(data_frame[variable_name])
    degree = int(degree)
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

    ## TODO I OVERWROTE CHEBYS WITH REGULARS
    #for i in np.arange(0,degree):
    #    X[:,i] = np.power(x,i+1)



    # generate names for each column
    col_names = []
    for i in np.arange(0,degree):
        col_names.append(variable_name + "^" + str(i+1))
    return X, col_names





















